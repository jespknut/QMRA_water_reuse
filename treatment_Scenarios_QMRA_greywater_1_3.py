### QMRA Simulator — CSV-Driven Version

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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
    'Rooftop rainwater': 0.01
}

# --- Helper Functions ---

EXPOSURE_ROUTES = {
    'ingestion': lambda v, f: v * f,
    'inhalation': lambda v, f: v * f * 0.1,
    'dermal': lambda v, f: v * f * 0.01
}

def get_raw_concentration(pathogen_name, matrix='Greywater'):
    row = raw_conc_df[
        (raw_conc_df['Pathogen'] == pathogen_name) &
        (raw_conc_df['Matrix'] == matrix)
    ]
    if not row.empty:
        return float(row.iloc[0]['Raw_Concentration'])
    else:
        return None

def sample_log_removal_sequence(treatments, pathogen_type, raw_concentration):
    total_removal = 0
    current_conc = raw_concentration
    for name, _, status in treatments:
        profiles = TREATMENT_EFFECTIVENESS_DIST.get(name, {}).get(pathogen_type, [])
        selected = next((p for p in profiles if current_conc <= p['max_conc']), profiles[-1]) if profiles else {'mean': 0, 'sd': 0}
        sampled = np.random.normal(selected['mean'], selected['sd'])
        effective_removal = sampled if status == 'normal' else sampled * 0.25 if status == 'degraded' else 0
        total_removal += effective_removal
        current_conc = max(current_conc / (10 ** effective_removal), 1e-3)
    return total_removal

def compute_probability_of_infection(dose, pathogen):
    if pathogen['Model_Type'] == 'Beta-Poisson':
        return 1 - (1 + dose / pathogen['Parameter_N50']) ** (-pathogen['Parameter_alpha'])
    else:
        k = np.random.normal(pathogen['Parameter_k'], 0.1 * pathogen['Parameter_k'])
        return 1 - np.exp(-k * dose)

def calculate_dose(conc, vol, pathogen):
    route = pathogen.get('Exposure_route', 'ingestion')
    absorption = 1.0
    route_func = EXPOSURE_ROUTES.get(route, EXPOSURE_ROUTES['ingestion'])
    return conc * route_func(vol, absorption)

def simulate_mc(activity, pathogen, treatment_steps=None, scenario_name="", source='Greywater', n_iter=10000, population=200):
    pathogen_type = pathogen.get('Type', 'bacteria')
    raw_conc = get_raw_concentration(pathogen['Pathogen'], matrix=source)
    if raw_conc is None:
        raise ValueError(f"No raw concentration data for pathogen {pathogen['Pathogen']} in matrix {source}")

    raw_conc *= SOURCE_WATER.get(source, 1.0)
    log_removal = sample_log_removal_sequence(treatment_steps, pathogen_type, raw_conc) if treatment_steps else 0
    adjusted_conc = max(raw_conc / (10 ** log_removal), 1e-6)
    mean_conc = np.log(adjusted_conc)

    conc = np.random.lognormal(mean=mean_conc, sigma=0.5, size=n_iter)
    vol = np.random.lognormal(mean=np.log(activity['Volume_L']), sigma=0.5, size=n_iter)
    vol = np.clip(vol, a_min=1e-8, a_max=None)

    freq = np.random.normal(loc=activity['Frequency_per_year'], scale=activity['Frequency_per_year'] * 0.1, size=n_iter)
    freq = np.clip(freq, a_min=1, a_max=None).astype(int)

    dose = calculate_dose(conc, vol, pathogen)
    p_inf = np.array([compute_probability_of_infection(d, pathogen) for d in dose])
    annual_risk = 1 - np.power((1 - p_inf), freq)

    expected_cases = annual_risk * population
    symptomatic_cases = expected_cases * pathogen['Symptomatic_Percent'] / 100
    hospital_days = expected_cases * pathogen['Hospitalized_Percent'] / 100 * pathogen['Hospital_days']
    daly = expected_cases * pathogen['DALY_per_case']

    return pd.DataFrame({
        'Pathogen': pathogen['Pathogen'],
        'Activity': activity['Activity'],
        'Scenario': scenario_name,
        'Source': source,
        'Annual risk': annual_risk,
        'Expected cases': expected_cases,
        'Sick days': symptomatic_cases * pathogen['Sick_days'],
        'Hospital days': hospital_days,
        'DALY': daly,
        'Dose': dose,
        'P_inf': p_inf
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
    print("\nSummary by scenario:")
    print(all_results.groupby("Scenario")[['Expected cases', 'Sick days', 'Hospital days', 'DALY']].mean())
    plot_scenario_matrix(all_results, 'Expected cases', f"Expected cases: {pathogen_name} at {activity_name}")
    plot_scenario_matrix(all_results, 'P_inf', f"Infection probability: {pathogen_name} at {activity_name}")
    plot_scenario_matrix(all_results, 'Annual risk', f"Annual risk: {pathogen_name} at {activity_name}")
    return all_results

def run_total_risk(activity_name, treatment_steps, scenario_name="All working", source='Greywater', n_iter=10000):
    activity = next(a for a in ACTIVITIES if a['Activity'] == activity_name)
    all_daly = []
    all_cases = []
    labels = []

    for pathogen in PATHOGENS:
        df = simulate_mc(activity, pathogen, treatment_steps=treatment_steps, scenario_name=scenario_name, source=source, n_iter=n_iter)
        all_daly.append(df['DALY'].reset_index(drop=True))
        all_cases.append(df['Expected cases'].reset_index(drop=True))
        labels.append(pathogen['Pathogen'])

    total_daly = pd.concat(all_daly, axis=1)
    total_daly.columns = labels
    total_cases = pd.concat(all_cases, axis=1)
    total_cases.columns = labels

    summary = pd.DataFrame({
        'Total DALY': total_daly.sum(axis=1),
        'Total cases': total_cases.sum(axis=1)
    })
    
    # Histogram of total DALY
    plt.figure(figsize=(8, 5))
    bins = np.logspace(np.log10(1e-9), np.log10(summary['Total DALY'].max() + 1e-8), 50)
    plt.hist(summary['Total DALY'], bins=bins, color='orchid', edgecolor='k', alpha=0.7)
    plt.xscale('log')

    p5 = np.percentile(summary['Total DALY'], 5)
    p50 = np.percentile(summary['Total DALY'], 50)
    p95 = np.percentile(summary['Total DALY'], 95)

    for p, label in zip([p5, p50, p95], ['5th percentile', 'Median', '95th percentile']):
        plt.axvline(p, linestyle='--', label=label)

    plt.axvspan(0, 1e-6, color='red', alpha=0.2, label='WHO threshold (1e-6)')
    plt.title(f"Total DALY distribution: {activity_name} ({scenario_name})")
    plt.xlabel("DALY/person-year (log scale)")
    plt.ylabel("Frequency")
    plt.legend()
    plt.tight_layout()
    plt.show()

    mean_dalys = total_daly.mean().sort_values(ascending=False)
    plt.figure(figsize=(10, 5))
    ax = mean_dalys.plot(kind='bar', color='teal', edgecolor='black')
    ax.set_yscale('log')
    plt.ylabel("Mean DALY contribution (log scale)")
    plt.title(f"Mean pathogen DALY contributions ({activity_name})")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

    return summary, total_daly

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
    
    run_treatment_scenarios("Shower", "Legionella pneumophila", treatment_scenarios, source='Greywater')
    summary, total_daly = run_total_risk("Shower", treatment_scenarios["All working"], scenario_name="All working", source='Greywater')
    summary['Total DALY'].to_csv("total_daly_script.csv", index=False)
