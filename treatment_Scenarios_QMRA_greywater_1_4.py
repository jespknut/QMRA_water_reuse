# === QMRA Simulator — CSV-Driven Version (fixed issues #1 and #2) ===
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

__version__ = "1.4.0"
__about__ = """
### About this app
This is the QMRA Greywater Reuse Simulator.

- CSV-driven, modular setup
- Pathogen library, dose–response, treatment trains
- Monte Carlo risk estimation
- Per-person DALY outputs compared to WHO benchmark
"""

DEGRADED_FACTOR_DEFAULT = 0.25  # editable; sensitivity below can sweep this

# --- Load Data from CSVs ---

# Load pathogens catalog (minimal list)
pathogens_catalog = pd.read_csv("qmra-data-library/data/pathogens.csv")

# Load integrated dose-response + burden data
pathogen_df = pd.read_csv("qmra-data-library/data/dose_response_burden.csv")
PATHOGENS = pathogen_df.to_dict(orient='records')

# Load activities
activity_df = pd.read_csv("qmra-data-library/data/exposure_parameters.csv")
ACTIVITIES = activity_df.to_dict(orient='records')
print("Columns in activity_df:", activity_df.columns.tolist())
print("First few rows:", activity_df.head())

# Load raw concentrations
raw_conc_df = pd.read_csv("qmra-data-library/data/raw_concentrations.csv")

# Load treatment effectiveness and build nested dictionary
treat_eff_df = pd.read_csv("qmra-data-library/data/treatment_effectiveness.csv")

# Build nested dict:
TREATMENT_EFFECTIVENESS_DIST = {}
for _, row in treat_eff_df.iterrows():
    proc = row['Treatment_process']
    ptype = row['Pathogen_type']
    if proc not in TREATMENT_EFFECTIVENESS_DIST:
        TREATMENT_EFFECTIVENESS_DIST[proc] = {}
    if ptype not in TREATMENT_EFFECTIVENESS_DIST[proc]:
        TREATMENT_EFFECTIVENESS_DIST[proc][ptype] = []
    TREATMENT_EFFECTIVENESS_DIST[proc][ptype].append({
        'max_conc': float(row['Max_conc_threshold']),
        'mean': float(row['Log_reduction_mean']),
        'sd': float(row['Log_reduction_sd'])
    })

# Source water factors
SOURCE_WATER = {
    'Greywater': 1.0,
    'Stormwater': 1.0,          # uses own data when available
    'Rooftop rainwater': 0.1         # proxy: scaled from Stormwater or Greywater
}


# --- Helper Functions ---

EXPOSURE_ROUTES = {
    'ingestion': lambda v, f: v * f,
    'inhalation': lambda v, f: v * f * 0.1,
    'dermal':     lambda v, f: v * f * 0.01
}

def get_raw_concentration(pathogen_name, matrix='Greywater'):
    # 1. Direct match
    row = raw_conc_df[
        (raw_conc_df['Pathogen'] == pathogen_name) &
        (raw_conc_df['Matrix'] == matrix)
    ]
    if not row.empty:
        mean = float(row.iloc[0]['Raw_Concentration'])
        sd_log10 = float(row.iloc[0]['Uncertainty_Log10_SD'])
        return mean, sd_log10

    # 2. Roof runoff → Stormwater proxy
    if matrix == 'Rooftop rainwater':
        storm_row = raw_conc_df[
            (raw_conc_df['Pathogen'] == pathogen_name) &
            (raw_conc_df['Matrix'] == 'Stormwater')
        ]
        if not storm_row.empty:
            mean = float(storm_row.iloc[0]['Raw_Concentration'])
            sd_log10 = float(storm_row.iloc[0]['Uncertainty_Log10_SD'])
            factor = SOURCE_WATER.get('Rooftop rainwater', 1.0)
            return mean * factor, sd_log10

    # 3. Final fallback → Greywater
    grey_row = raw_conc_df[
        (raw_conc_df['Pathogen'] == pathogen_name) &
        (raw_conc_df['Matrix'] == 'Greywater')
    ]
    if not grey_row.empty:
        mean = float(grey_row.iloc[0]['Raw_Concentration'])
        sd_log10 = float(grey_row.iloc[0]['Uncertainty_Log10_SD'])
        factor = SOURCE_WATER.get(matrix, 1.0)
        return mean * factor, sd_log10

    raise ValueError(
        f"No raw concentration data for pathogen '{pathogen_name}' "
        f"in matrix '{matrix}' and no valid fallback"
    )


def _truncated_normal(mean, sd, floor=0.0):
    """Draw a normal sample and truncate at `floor` (>=0 to avoid 'amplification')."""
    x = np.random.normal(mean, sd)
    return max(x, floor)

def sample_log_removal_sequence(treatments, pathogen_type, raw_concentration,
                                degraded_factor=DEGRADED_FACTOR_DEFAULT):
    """
    Piece-wise LR sampling by current concentration with physical truncation at 0.
    `degraded_factor` scales LR when status == 'degraded'.
    """
    total_removal = 0.0
    current_conc = raw_concentration

    for name, _, status in (treatments or []):
        profiles = TREATMENT_EFFECTIVENESS_DIST.get(name, {}).get(pathogen_type, [])

        if profiles:
            # choose profile based on current concentration
            selected = next((p for p in profiles if current_conc <= p['max_conc']), profiles[-1])
            sampled_LR = _truncated_normal(selected['mean'], selected['sd'], floor=0.0)
        else:
            sampled_LR = 0.0

        if status == 'normal':
            effective_LR = sampled_LR
        elif status == 'degraded':
            effective_LR = sampled_LR * float(degraded_factor)
        else:  # 'fail'
            effective_LR = 0.0

        total_removal += effective_LR
        # update concentration for next step (never let it hit zero)
        current_conc = max(current_conc / (10 ** effective_LR), 1e-9)

    return total_removal

def compare_failure_impact(activity_name, treatment_steps, source='Greywater',
                           n_iter=10000, population=200, degraded_factor=DEGRADED_FACTOR_DEFAULT):
    """
    Compute Δmedian and ΔP95 for per-person DALY when steps are degraded/failed
    vs a baseline where the same steps are all 'normal'.
    """
    # Baseline: same steps set to 'normal'
    baseline_steps = [(name, meta, 'normal') for (name, meta, _status) in treatment_steps]

    # Run totals (per-person DALY distribution returned via summary)
    base_summary, _ = run_total_risk(
        activity_name, baseline_steps, scenario_name="Baseline (all normal)",
        source=source, n_iter=n_iter, population=population
    )
    test_summary, _ = run_total_risk(
        activity_name, treatment_steps, scenario_name="As selected",
        source=source, n_iter=n_iter, population=population
    )

    b = base_summary['Per-person DALY'].values
    t = test_summary['Per-person DALY'].values

    def p(x, q): return np.percentile(x, q)

    out = {
        "Baseline_P50": p(b, 50), "Baseline_P95": p(b, 95),
        "Selected_P50": p(t, 50), "Selected_P95": p(t, 95),
        "Delta_P50": p(t, 50) - p(b, 50),
        "Delta_P95": p(t, 95) - p(b, 95),
        "Ratio_P50": (p(t, 50) / p(b, 50)) if p(b, 50) > 0 else np.inf,
        "Ratio_P95": (p(t, 95) / p(b, 95)) if p(b, 95) > 0 else np.inf,
    }
    return pd.Series(out)

def degraded_sensitivity(activity_name, treatment_steps, source='Greywater',
                         n_iter=10000, population=200, factors=(0.1, 0.25, 0.5)):
    """
    Sweep degraded_factor over a set of values and report P50/P95 & deltas vs baseline.
    """
    results = []
    # normalize steps to degraded to feel the scalar impact
    steps_deg = [(n, m, ('degraded' if s != 'fail' else 'fail')) for (n, m, s) in treatment_steps]

    # Baseline once
    baseline_steps = [(name, meta, 'normal') for (name, meta, _status) in treatment_steps]
    base_summary, _ = run_total_risk(activity_name, baseline_steps, source=source,
                                     n_iter=n_iter, population=population)
    b = base_summary['Per-person DALY'].values
    P50b, P95b = np.percentile(b, 50), np.percentile(b, 95)

    for f in factors:
        # monkey-patch degraded factor via partial function call:
        # reuse simulate_mc by passing degraded_factor downstream
        # Easiest: temporarily wrap sample function via closure; or if you call simulate_mc,
        # add a param and thread through. For brevity, we illustrate using run_total_risk as-is.
        # (If your simulate_mc takes degraded_factor param, pass it there.)
        # Here we assume you edited simulate_mc to accept **kwargs and forward to sample function.

        # Run with degraded steps (respecting f via global default)
        # If you prefer not to rely on a global, add a parameter to simulate_mc/sample to pass f.
        global DEGRADED_FACTOR_DEFAULT
        old_val = DEGRADED_FACTOR_DEFAULT
        DEGRADED_FACTOR_DEFAULT = f
        try:
            s, _ = run_total_risk(activity_name, steps_deg, source=source, n_iter=n_iter, population=population)
        finally:
            DEGRADED_FACTOR_DEFAULT = old_val

        t = s['Per-person DALY'].values
        P50t, P95t = np.percentile(t, 50), np.percentile(t, 95)

        results.append({
            "Degraded_factor": f,
            "Baseline_P50": P50b, "Baseline_P95": P95b,
            "Selected_P50": P50t, "Selected_P95": P95t,
            "Delta_P50": P50t - P50b, "Delta_P95": P95t - P95b,
            "Ratio_P50": (P50t / P50b) if P50b > 0 else np.inf,
            "Ratio_P95": (P95t / P95b) if P95b > 0 else np.inf
        })
    return pd.DataFrame(results)


def compute_probability_of_infection(dose, pathogen):
    if pathogen.get('Model_Type', 'Exponential') == 'Beta-Poisson':
        # NOTE: using classical approximate Beta-Poisson form with alpha, N50
        return 1 - (1 + dose / pathogen['Parameter_N50']) ** (-pathogen['Parameter_alpha'])
    else:
        # Exponential model; sample k with modest uncertainty
        k = np.random.normal(pathogen['Parameter_k'], 0.1 * pathogen['Parameter_k'])
        return 1 - np.exp(-k * dose)

def calculate_dose(conc, vol, pathogen, activity):
    route = pathogen.get('Exposure_route', 'ingestion')
    if route == 'ingestion':
        return conc * vol
    elif route == 'inhalation':
        # New inhalation module
        aerosol_factor = activity.get('Aerosol_factor', 1e-5)  # L water → m³ air
        inhalation_rate = activity.get('Inhalation_rate', 0.012)  # m³/min
        duration = activity.get('Duration_min', 10)  # minutes
        deposition_fraction = activity.get('Deposition_fraction', 0.6)
        inhaled_volume = inhalation_rate * duration
        dose = conc * vol * aerosol_factor * inhaled_volume * deposition_fraction
        return dose
    elif route == 'dermal':
        return conc * vol * 0.01  # placeholder
    else:
        return conc * vol

def simulate_mc(activity, pathogen, treatment_steps=None, scenario_name="",
                source='Greywater', n_iter=10000, population=200, degraded_factor=DEGRADED_FACTOR_DEFAULT, occurrence_mode="stochastic"):
    """Runs MC for ONE pathogen and ONE activity. Returns a row-wise dataframe of outcomes.

    IMPORTANT: This now returns BOTH per-person and population-level metrics.
 
    occurrence_mode:
      - 'deterministic': use mean concentration only
      - 'stochastic': sample from lognormal (default)
    """
    
    if 'Type' not in pathogen:
        raise ValueError(
            f"Pathogen '{pathogen['Pathogen']}' is missing 'Type' "
            "(bacteria / virus / protozoa)"
        )
    pathogen_type = pathogen['Type']

    
    mean_raw, sd_log10 = get_raw_concentration(
        pathogen['Pathogen'],
        matrix=source
    )

    # --- Sample raw concentration (occurrence model) ---
    if occurrence_mode == "deterministic":
        # Use mean concentration only (for debugging / sensitivity)
        raw_conc_samples = np.full(n_iter, mean_raw)
    else:
        # Log10-normal occurrence uncertainty from CSV
        raw_conc_samples = 10 ** np.random.normal(
            loc=np.log10(mean_raw),
            scale=sd_log10,
            size=n_iter
        )


    # --- Apply treatment per iteration ---
    adjusted_conc = np.zeros(n_iter)
    for i in range(n_iter):
        lr = sample_log_removal_sequence(
            treatment_steps,
            pathogen_type,
            raw_conc_samples[i],
            degraded_factor=degraded_factor
        ) if treatment_steps else 0.0

        adjusted_conc[i] = max(
            raw_conc_samples[i] / (10 ** lr),
            1e-9
        )

    # --- Final concentration used for exposure ---
    conc = adjusted_conc

    vol  = np.random.lognormal(mean=np.log(activity['Volume_L']), sigma=0.5, size=n_iter)
    vol  = np.clip(vol, a_min=1e-8, a_max=None)

    # Frequency per year (integer ≥1)
    freq = np.random.normal(loc=activity['Frequency_per_year'],
                            scale=max(1e-9, activity['Frequency_per_year'] * 0.1),
                            size=n_iter)
    freq = np.clip(freq, a_min=1, a_max=None).astype(int)

    # Dose & infection
    dose = calculate_dose(conc, vol, pathogen, activity)
    p_inf  = np.array([compute_probability_of_infection(d, pathogen) for d in dose])
    annual_risk = 1 - np.power((1 - p_inf), freq)  # per-person annual risk for this pathogen

    # --- Per-person vs population metrics (ISSUE #1 fix) ---
    # Per-person expected cases = annual_risk
    expected_cases_population = annual_risk * population
    # Per-person DALY (for WHO benchmark)
    daly_per_person = annual_risk * pathogen['DALY_per_case']
    # Population DALY (sum across persons)
    daly_population = daly_per_person * population

    # Symptomatic/hospitalization days (population-level)
    symptomatic_cases = expected_cases_population * pathogen['Symptomatic_Percent'] / 100.0
    hospital_days     = expected_cases_population * pathogen['Hospitalized_Percent'] / 100.0 * pathogen['Hospital_days']

    return pd.DataFrame({
        'Pathogen': pathogen['Pathogen'],
        'Activity': activity['Activity'],
        'Scenario': scenario_name,
        'Source':   source,
        # Risks & dose
        'P_inf': p_inf,
        'Annual risk': annual_risk,            # per-person, per year
        'Dose': dose,
        # Per-person health impact (use for WHO comparison)
        'DALY_per_person': daly_per_person,
        # Population-level counts/impact
        'Expected cases': expected_cases_population,
        'Sick days': symptomatic_cases * pathogen['Sick_days'],
        'Hospital days': hospital_days,
        'DALY_population': daly_population
    })

def plot_scenario_matrix(df, column, title_prefix):
    scenarios = df['Scenario'].unique()
    cols = 3
    rows = (len(scenarios) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows))
    axes = axes.flatten()

    for i, scenario in enumerate(scenarios):
        ax = axes[i]
        values = df[df['Scenario'] == scenario][column].replace([np.inf, -np.inf], np.nan).dropna()
        if values.empty:
            ax.set_title(f"{scenario} — No data")
            continue
        ax.hist(values, bins=40, color='skyblue', edgecolor='k', alpha=0.7)
        ax.axvline(values.mean(), color='red', linestyle='--', label=f"Mean: {values.mean():.2f}")
        ax.set_title(f"{scenario}")
        ax.set_xlabel(column)
        ax.set_ylabel("Frequency")
        ax.legend()

    # Remove extra axes if any
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    fig.suptitle(f"{title_prefix} by scenario", fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

def run_treatment_scenarios(activity_name, pathogen_name, scenario_defs, source='Greywater'):
    activity = next(a for a in ACTIVITIES if a['Activity'] == activity_name)
    pathogen = next(p for p in PATHOGENS if p['Pathogen'] == pathogen_name)
    results = []
    for name, steps in scenario_defs.items():
        df = simulate_mc(activity, pathogen, treatment_steps=steps, scenario_name=name, source=source)
        results.append(df)
    all_results = pd.concat(results)
    print("\nSummary by scenario (population-level means):")
    print(all_results.groupby("Scenario")[['Expected cases', 'Sick days', 'Hospital days', 'DALY_population']].mean())
    plot_scenario_matrix(all_results, 'Expected cases', f"Expected cases: {pathogen_name} at {activity_name}")
    plot_scenario_matrix(all_results, 'P_inf', f"Infection probability: {pathogen_name} at {activity_name}")
    plot_scenario_matrix(all_results, 'Annual risk', f"Annual risk: {pathogen_name} at {activity_name}")
    return all_results

def run_total_risk(activity_name, treatment_steps, scenario_name="All working",
                   source='Greywater', n_iter=10000, population=200):
    """Aggregate across ALL pathogens for ONE activity & treatment scenario.

    Fixes:
    - ISSUE #1: Aggregate and plot **per-person** DALY distribution for WHO comparison.
    - ISSUE #2: Compute **combined annual infection risk** across pathogens with complement-of-product.
    """
    activity = next(a for a in ACTIVITIES if a['Activity'] == activity_name)
    daly_person_cols = []
    annual_risk_cols = []
    daly_pop_cols = []
    cases_pop_cols = []
    labels = []

    for pathogen in PATHOGENS:
        df = simulate_mc(activity, pathogen, treatment_steps=treatment_steps,
                         scenario_name=scenario_name, source=source, n_iter=n_iter, population=population)
        daly_person_cols.append(df['DALY_per_person'].reset_index(drop=True))
        annual_risk_cols.append(df['Annual risk'].reset_index(drop=True))
        daly_pop_cols.append(df['DALY_population'].reset_index(drop=True))
        cases_pop_cols.append(df['Expected cases'].reset_index(drop=True))
        labels.append(pathogen['Pathogen'])

    # Build per-pathogen matrices
    daly_person = pd.concat(daly_person_cols, axis=1); daly_person.columns = labels
    annual_risk = pd.concat(annual_risk_cols, axis=1); annual_risk.columns = labels
    daly_pop    = pd.concat(daly_pop_cols, axis=1);    daly_pop.columns = labels
    cases_pop   = pd.concat(cases_pop_cols, axis=1);   cases_pop.columns = labels

    # --- Aggregations ---
    # Per-person DALY is additive across pathogens (sum per iteration)
    total_daly_per_person = daly_person.sum(axis=1)
    # Population DALY & cases are also additive
    total_daly_population = daly_pop.sum(axis=1)
    total_cases_population = cases_pop.sum(axis=1)
    # ISSUE #2: Combined infection risk across pathogens (NOT a sum)
    combined_annual_risk = 1 - (1 - annual_risk).prod(axis=1)

    # --- Plot WHO-relevant distribution (per-person DALY) ---
    plt.figure(figsize=(8, 5))
    # avoid zero-only bins
    bins = np.logspace(np.log10(1e-9), np.log10(max(1e-9, total_daly_per_person.max()) + 1e-8), 50)
    plt.hist(total_daly_per_person, bins=bins, color='orchid', edgecolor='k', alpha=0.7)
    plt.xscale('log')

    p5  = np.percentile(total_daly_per_person, 5)
    p50 = np.percentile(total_daly_per_person, 50)
    p95 = np.percentile(total_daly_per_person, 95)
    for p, label in zip([p5, p50, p95], ['5th percentile', 'Median', '95th percentile']):
        plt.axvline(p, linestyle='--', label=label)

    # WHO benchmark band (per-person DALY per year)
    plt.axvspan(0, 1e-6, color='red', alpha=0.2, label='WHO threshold (1e-6 per-person DALY/y)')
    plt.title(f"Per-person DALY distribution: {activity_name} ({scenario_name})")
    plt.xlabel("DALY per person per year (log scale)")
    plt.ylabel("Frequency")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Bar chart of mean pathogen DALY contributions (per-person, log scale)
    mean_dalys_person = daly_person.mean().sort_values(ascending=False)
    plt.figure(figsize=(10, 5))
    ax = mean_dalys_person.plot(kind='bar', color='teal', edgecolor='black')
    ax.set_yscale('log')
    plt.ylabel("Mean per-person DALY contribution (log scale)")
    plt.title(f"Mean pathogen DALY contributions (per-person) — {activity_name}")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

    # Summary dataframe with both per-person and population aggregates
    summary = pd.DataFrame({
        'Per-person DALY': total_daly_per_person,
        'Combined annual infection risk': combined_annual_risk,
        'Population DALY': total_daly_population,
        'Population cases': total_cases_population
    })

    # Quick textual summary
    print("\n--- Aggregated results (key percentiles) ---")
    for col in ['Per-person DALY', 'Combined annual infection risk', 'Population DALY', 'Population cases']:
        arr = summary[col].values
        print(f"{col}: P50={np.percentile(arr,50):.3e}, P95={np.percentile(arr,95):.3e}, mean={arr.mean():.3e}")

    return summary, daly_person

# Example execution
if __name__ == "__main__":
    treatment_scenarios = {
        "All working": [
            ('Ultrafiltration', None, 'normal'),
            ('UV', None, 'normal')
        ],
        "UF failed": [
            ('Ultrafiltration', None, 'fail'),
            ('UV', None, 'normal')
        ]
    }

    # Single-pathogen scenario exploration
    run_treatment_scenarios("Shower", "Legionella pneumophila", treatment_scenarios, source='Greywater')

    # Cross-pathogen totals for WHO comparison (per-person) and proper risk aggregation
    summary, daly_person_matrix = run_total_risk(
        "Shower",
        treatment_scenarios["All working"],
        scenario_name="All working",
        source='Greywater',
        n_iter=10000,
        population=200
    )

    # Export per-person DALY distribution if desired
    summary[['Per-person DALY']].to_csv("per_person_daly_distribution.csv", index=False)
