### QMRA Simulator — Finalized with Concentration-Aware Treatment and Beta-Poisson

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# --- Data ---

SOURCE_WATER = {
    'Greywater': 1.0,
    'Rooftop rainwater': 0.01  # Assumes 100x lower baseline pathogen levels
}

PATHOGENS = [
    {'Name': 'Escherichia coli', 'Raw concentration': 70000, 'alpha': 0.155, 'N50': 2.11e6, 'Type': 'bacteria', 'Symptomatic (%)': 32, 'Sick days': 3, 'Hospitalized (%)': 1, 'Hospital days': 2, 'DALY per case': 1e-4},
    {'Name': 'Giardia', 'Raw concentration': 1.2, 'alpha': 0.359, 'N50': 50, 'Type': 'protozoa', 'Symptomatic (%)': 80, 'Sick days': 7, 'Hospitalized (%)': 2, 'Hospital days': 3, 'DALY per case': 4e-4},
    {'Name': 'Salmonella spp.', 'Raw concentration': 574000, 'k': 3.97e-06, 'Type': 'bacteria', 'Symptomatic (%)': 80, 'Sick days': 4, 'Hospitalized (%)': 2, 'Hospital days': 2, 'DALY per case': 2e-4},
    {'Name': 'Shigella spp.', 'Raw concentration': 50000, 'k': 0.0315, 'Type': 'bacteria', 'Symptomatic (%)': 85, 'Sick days': 5, 'Hospitalized (%)': 3, 'Hospital days': 2, 'DALY per case': 3e-4},
    {'Name': 'Rotavirus', 'Raw concentration': 1e7, 'alpha': 0.253, 'N50': 6.17e5, 'Type': 'virus', 'Symptomatic (%)': 95, 'Sick days': 6, 'Hospitalized (%)': 2, 'Hospital days': 2, 'DALY per case': 4e-4},
    {'Name': 'Cryptosporidium parvum', 'Raw concentration': 11000, 'alpha': 1.14e-01, 'N50': 4.55e2, 'Type': 'protozoa', 'Symptomatic (%)': 70, 'Sick days': 7, 'Hospitalized (%)': 1, 'Hospital days': 3, 'DALY per case': 5e-4},
    {'Name': 'Mycobacterium avium', 'Raw concentration': 2800, 'k': 6.93e-04, 'Type': 'bacteria', 'Symptomatic (%)': 5, 'Sick days': 10, 'Hospitalized (%)': 1, 'Hospital days': 5, 'DALY per case': 2e-3},
    {'Name': 'Pseudomonas aeruginosa', 'Raw concentration': 17783.9, 'k': 1.05e-4, 'Type': 'bacteria', 'Symptomatic (%)': 10, 'Sick days': 4, 'Hospitalized (%)': 1, 'Hospital days': 2, 'DALY per case': 1e-3},
    {'Name': 'Adenovirus', 'Raw concentration': 1e6, 'alpha': 0.4172, 'N50': 1.12e5, 'Type': 'virus', 'Symptomatic (%)': 75, 'Sick days': 5, 'Hospitalized (%)': 2, 'Hospital days': 2, 'DALY per case': 2e-4},
    {'Name': 'Enterovirus', 'Raw concentration': 1e6, 'k': 3.74e-03, 'Type': 'virus', 'Symptomatic (%)': 70, 'Sick days': 4, 'Hospitalized (%)': 1, 'Hospital days': 1, 'DALY per case': 1e-4},
    {'Name': 'Hepatit A-virus', 'Raw concentration': 1e6, 'alpha': 0.253, 'N50': 1e6, 'Type': 'virus', 'Symptomatic (%)': 90, 'Sick days': 14, 'Hospitalized (%)': 10, 'Hospital days': 5, 'DALY per case': 8e-4},
    {'Name': 'Legionella pneumophila', 'Raw concentration': 4500, 'k': 5.99e-02, 'Type': 'bacteria', 'Symptomatic (%)': 5, 'Sick days': 10, 'Hospitalized (%)': 2, 'Hospital days': 10, 'DALY per case': 1e-3, 'Exposure route': 'inhalation'}
]
    

ACTIVITIES = [
    {'Name': 'Shower', 'Volume (L)': 0.001, 'Frequency/year': 365},
    {'Name': 'Toilet flushing', 'Volume (L)': 0.00001, 'Frequency/year': 1000},
    {'Name': 'Handwashing', 'Volume (L)': 0.00001, 'Frequency/year': 400},
    {'Name': 'Irrigation', 'Volume (L)': 0.0005, 'Frequency/year': 60},
    {'Name': 'Dishwashing', 'Volume (L)': 0.0005, 'Frequency/year': 300},
    {'Name': 'Spray irrigation', 'Volume (L)': 0.0002, 'Frequency/year': 240}  # Inhalation or dermal possible
]

TREATMENT_EFFECTIVENESS_DIST = {
    'Coarse filter': {
        'virus': [{'max_conc': 1e8, 'mean': 0.5, 'sd': 0.1}],
        'bacteria': [{'max_conc': 1e8, 'mean': 0.5, 'sd': 0.1}],
        'protozoa': [{'max_conc': 1e8, 'mean': 0.2, 'sd': 0.05}]
    },
    'Ultrafiltration': {
        'virus': [{'max_conc': 1e8, 'mean': 3, 'sd': 0.3}],
        'bacteria': [{'max_conc': 1e8, 'mean': 4, 'sd': 0.3}],
        'protozoa': [{'max_conc': 1e8, 'mean': 5, 'sd': 0.3}]
    },
    'Nanofiltration': {
        'virus': [{'max_conc': 1e4, 'mean': 3, 'sd': 0.3}, {'max_conc': 1e8, 'mean': 5, 'sd': 0.3}],
        'bacteria': [{'max_conc': 1e4, 'mean': 4, 'sd': 0.3}, {'max_conc': 1e8, 'mean': 6, 'sd': 0.3}],
        'protozoa': [{'max_conc': 1e4, 'mean': 4, 'sd': 0.2}, {'max_conc': 1e8, 'mean': 5, 'sd': 0.3}]
    },
    'UV': {
        'virus': [{'max_conc': 1e8, 'mean': 4, 'sd': 0.2}],
        'bacteria': [{'max_conc': 1e8, 'mean': 3, 'sd': 0.2}],
        'protozoa': [{'max_conc': 1e8, 'mean': 1.5, 'sd': 0.1}]
    },
    'Chlorination': {
        'virus': [{'max_conc': 1e4, 'mean': 3, 'sd': 0.3}, {'max_conc': 1e8, 'mean': 5, 'sd': 0.3}],
        'bacteria': [{'max_conc': 1e4, 'mean': 2, 'sd': 0.2}, {'max_conc': 1e8, 'mean': 3, 'sd': 0.2}],
        'protozoa': [{'max_conc': 1e4, 'mean': 1, 'sd': 0.2}, {'max_conc': 1e8, 'mean': 1.5, 'sd': 0.2}]
    },
    'Sand filter': {
        'virus': [{'max_conc': 1e8, 'mean': 1.2, 'sd': 0.4}],
        'bacteria': [{'max_conc': 1e8, 'mean': 2.0, 'sd': 0.5}],
        'protozoa': [{'max_conc': 1e8, 'mean': 2.5, 'sd': 0.5}]
    },
    'GAC filter': {
        'virus': [{'max_conc': 1e8, 'mean': 1.5, 'sd': 0.5}],
        'bacteria': [{'max_conc': 1e8, 'mean': 2.5, 'sd': 0.6}],
        'protozoa': [{'max_conc': 1e8, 'mean': 2.0, 'sd': 0.5}]
    },
    'Raised plant bed': {
        'virus': [{'max_conc': 1e8, 'mean': 1.2, 'sd': 0.4}],
        'bacteria': [{'max_conc': 1e8, 'mean': 2.5, 'sd': 0.6}],
        'protozoa': [{'max_conc': 1e8, 'mean': 2.0, 'sd': 0.5}]
    }
}

TREATMENT_EFFECTIVENESS_DIST['Ceramic membrane UF'] = {
    'virus': [{'max_conc': 1e8, 'mean': 3, 'sd': 0.3}],
    'bacteria': [{'max_conc': 1e8, 'mean': 4, 'sd': 0.3}],
    'protozoa': [{'max_conc': 1e8, 'mean': 5, 'sd': 0.2}]
}

# --- Helper Functions ---

EXPOSURE_ROUTES = {
    'ingestion': lambda v, f: v * f,
    'inhalation': lambda v, f: v * f * 0.1,
    'dermal': lambda v, f: v * f * 0.01
}

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
    if 'alpha' in pathogen and 'N50' in pathogen:
        return 1 - (1 + dose / pathogen['N50']) ** (-pathogen['alpha'])
    k = np.random.normal(pathogen['k'], 0.1 * pathogen['k'])
    return 1 - np.exp(-k * dose)

def calculate_dose(conc, vol, pathogen):
    route = pathogen.get('Exposure route', 'ingestion')
    absorption = pathogen.get('Absorption Fraction', 1.0)
    route_func = EXPOSURE_ROUTES.get(route, EXPOSURE_ROUTES['ingestion'])
    return conc * route_func(vol, absorption)

def simulate_mc(activity, pathogen, treatment_steps=None, scenario_name="", source='Greywater', n_iter=10000, population=200):
    pathogen_type = pathogen.get('Type', 'bacteria')
    raw_conc = pathogen['Raw concentration'] * SOURCE_WATER.get(source, 1.0)
    log_removal = sample_log_removal_sequence(treatment_steps, pathogen_type, raw_conc) if treatment_steps else 0
    adjusted_conc = max(raw_conc / (10 ** log_removal), 1e-6)
    mean_conc = np.log(adjusted_conc)

    conc = np.random.lognormal(mean=mean_conc, sigma=0.5, size=n_iter)
    vol = np.random.lognormal(mean=np.log(activity['Volume (L)']), sigma=0.5, size=n_iter)
    vol = np.clip(vol, a_min=1e-8, a_max=None)

    freq = np.random.normal(loc=activity['Frequency/year'], scale=activity['Frequency/year'] * 0.1, size=n_iter)
    freq = np.clip(freq, a_min=1, a_max=None).astype(int)

    dose = calculate_dose(conc, vol, pathogen)
    p_inf = np.array([compute_probability_of_infection(d, pathogen) for d in dose])
    annual_risk = 1 - np.power((1 - p_inf), freq)

    expected_cases = annual_risk * population
    symptomatic_cases = expected_cases * pathogen['Symptomatic (%)'] / 100
    hospital_days = expected_cases * pathogen['Hospitalized (%)'] / 100 * pathogen['Hospital days']
    daly = expected_cases * pathogen['DALY per case']

    return pd.DataFrame({
        'Pathogen': pathogen['Name'],
        'Activity': activity['Name'],
        'Scenario': scenario_name,
        'Source': source,
        'Annual risk': annual_risk,
        'Expected cases': expected_cases,
        'Sick days': symptomatic_cases * pathogen['Sick days'],
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
    activity = next(a for a in ACTIVITIES if a['Name'] == activity_name)
    pathogen = next(p for p in PATHOGENS if p['Name'] == pathogen_name)
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
    activity = next(a for a in ACTIVITIES if a['Name'] == activity_name)
    all_daly = []
    all_cases = []
    labels = []

    for pathogen in PATHOGENS:
        df = simulate_mc(activity, pathogen, treatment_steps=treatment_steps, scenario_name=scenario_name, source=source, n_iter=n_iter)
        all_daly.append(df['DALY'].reset_index(drop=True))
        all_cases.append(df['Expected cases'].reset_index(drop=True))
        labels.append(pathogen['Name'])

    total_daly = pd.concat(all_daly, axis=1)
    total_daly.columns = labels
    total_cases = pd.concat(all_cases, axis=1)
    total_cases.columns = labels

    summary = pd.DataFrame({
        'Total DALY': total_daly.sum(axis=1),
        'Total cases': total_cases.sum(axis=1)
    })
    
    # Histogram of total DALY with log-scale x-axis and percentile bands
    plt.figure(figsize=(8, 5))
    bins = np.logspace(np.log10(1e-9), np.log10(summary['Total DALY'].max() + 1e-8), 50)
    plt.hist(summary['Total DALY'], bins=bins,
             color='orchid', edgecolor='k', alpha=0.7)
    plt.xscale('log')

    # Percentiles
    p5 = np.percentile(summary['Total DALY'], 5)
    p50 = np.percentile(summary['Total DALY'], 50)
    p95 = np.percentile(summary['Total DALY'], 95)

    for p, label in zip([p5, p50, p95], ['5th percentile', 'Median', '95th percentile']):
        plt.axvline(p, linestyle='--', label=label)

    # WHO threshold as shaded region
    plt.axvspan(0, 1e-6, color='red', alpha=0.2, label='WHO threshold (1e-6)')

    plt.title(f"Total DALY distribution: {activity_name} ({scenario_name})")
    plt.xlabel("DALY/person-year (log scale)")
    plt.ylabel("Frequency")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Plot with log-scale x-axis for DALY contributions
    plt.figure(figsize=(10, 5))
    mean_dalys = total_daly.mean().sort_values(ascending=False)
    ax = mean_dalys.plot(kind='bar', color='teal', edgecolor='black')
    ax.set_yscale('log')
    plt.ylabel("Mean DALY contribution (log scale)")
    plt.title(f"Mean pathogen DALY contributions ({activity_name})")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()


    return summary, total_daly


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
    
