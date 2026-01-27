# 🧪 QMRA Greywater Reuse Simulator

## Overview
This tool is a **flexible, CSV-driven Quantitative Microbial Risk Assessment (QMRA) simulator** for **greywater reuse scenarios**.

It integrates:

- **Pathogen occurrence data**
- **Exposure parameters** (activity frequency, volume, route)
- **Dose–response models**
- **Treatment trains** with uncertainty, degraded, and failure states
- **Monte Carlo simulation** for probabilistic risk assessment

It outputs distributions of:

- Infection risk (per-person, per-year)
- Expected cases (population-level)
- DALYs (per-person and population)
- Symptomatic and hospitalized cases

Plots compare outcomes against the **WHO benchmark of 1×10⁻⁶ DALY/person-year**.

---

## 📂 File Structure

```
treatment_Scenarios_QMRA_greywater_1_3.py   # Main simulation script
qmra-data-library/
  ├── pathogens.csv                         # Pathogen metadata (name, type)
  ├── dose_response_burden.csv              # Dose-response params & DALY per case
  ├── exposure_parameters.csv               # Activity parameters (volume, freq, route)
  ├── raw_concentrations.csv                # Pathogen concentrations in raw greywater
  └── treatment_effectiveness.csv           # Treatment LRV distributions
```

---

## ⚙️ Core Functions

### `simulate_mc(activity, pathogen, treatment_steps, ...)`
Run a **Monte Carlo simulation** for one pathogen + one activity.  

**Inputs:**
- `activity`: dict from `exposure_parameters.csv`
- `pathogen`: dict from `dose_response_burden.csv`
- `treatment_steps`: list of `(process_name, _, status)`  
  - `status`: `"normal"`, `"degraded"`, `"fail"`
- `n_iter`: number of MC iterations
- `population`: population size

**Outputs:** `pandas.DataFrame` with columns:

| Column            | Description                                  |
|-------------------|----------------------------------------------|
| `P_inf`           | Per-exposure infection probability           |
| `Annual risk`     | Per-person annual infection risk             |
| `Expected cases`  | Population cases/year                        |
| `DALY_per_person` | DALY/person/year (for WHO comparison)        |
| `DALY_population` | Population DALY/year                         |
| `Sick days`       | Symptomatic case-days                        |
| `Hospital days`   | Hospitalization days                         |
| `Dose`            | Exposure dose (pathogen units)               |

---

### `run_treatment_scenarios(activity_name, pathogen_name, scenario_defs, source)`
Compare multiple **treatment scenarios** for a single pathogen + activity.

**Inputs:**
- `activity_name`: e.g., `"Shower"`
- `pathogen_name`: e.g., `"Legionella pneumophila"`
- `scenario_defs`: dict of scenarios, e.g.:

```python
treatment_scenarios = {
    "All working": [
        ("Ultrafiltration", None, "normal"),
        ("UV", None, "normal")
    ],
    "UF failed": [
        ("Ultrafiltration", None, "fail"),
        ("UV", None, "normal")
    ]
}
```

- `source`: `"Greywater"` or `"Rooftop rainwater"`

**Outputs:**
- Population-level summary table
- Histograms of expected cases, infection probability, annual risk
- Returns concatenated dataframe

---

### `run_total_risk(activity_name, treatment_steps, ...)`
Aggregate risk across **all pathogens** for one activity + treatment train.

✅ **Fixes implemented**:
- Infection risks combined using **complement-of-product**:  
  \[
  P_\text{total} = 1 - \prod_i (1 - P_i)
  \]
- DALYs and cases summed (additive quantities)

**Outputs:**
- Histogram of **per-person DALY** vs WHO benchmark (log scale)
- Bar chart of pathogen DALY contributions
- Percentile summary (P50, P95, mean) for:
  - Per-person DALY
  - Combined annual infection risk
  - Population DALY
  - Population cases
- Returns `(summary_df, daly_person_matrix)`

---

## ▶️ Example Usage

```python
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

    # Single-pathogen exploration
    run_treatment_scenarios("Shower", "Legionella pneumophila", treatment_scenarios)

    # Aggregated across all pathogens (WHO comparison)
    summary, daly_matrix = run_total_risk(
        "Shower",
        treatment_scenarios["All working"],
        scenario_name="All working",
        population=200
    )
```

---

## 📊 Key Assumptions

- Pathogen occurrence ~ **lognormal** (σ=0.5)
- Exposure volume ~ **lognormal** (σ=0.5)
- Frequency ~ **normal ±10%**
- Treatment LRVs ~ **normal distributions**, truncated ≥0
- **Degraded state** = 25% of nominal LRV
- **Failure state** = 0 LRV
- DALYs benchmarked as **per-person per year** against WHO 1×10⁻⁶

---

## 📈 Outputs

- **Histograms** of per-person DALY (log scale, WHO band shaded)
- **Scenario-wise distributions** of infection risk & cases
- **CSV export**: `per_person_daly_distribution.csv`

---

## 🔁 Reproducibility

- Monte Carlo sample size: **10,000** iterations (default)
- Random seeds can be fixed (`np.random.seed(...)`)
- All inputs configurable via **CSV libraries**

---

## 📝 Citation

If you use this tool in research, please cite:

> *QMRA Greywater Reuse Simulator (2025). A CSV-driven probabilistic tool for microbial risk assessment in decentralized water reuse.*
