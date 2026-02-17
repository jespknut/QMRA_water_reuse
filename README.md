# QMRA Water Reuse Simulator

A transparent, modular **Quantitative Microbial Risk Assessment (QMRA)** tool for evaluating health risks associated with **water reuse scenarios**, with explicit comparison to the WHO health-based target of **1×10⁻⁶ DALY per person-year**.

The tool is designed for **research, scenario analysis, and decision support**, with a strong emphasis on methodological clarity, reproducibility, and conservative risk aggregation.

---

## ✨ Key Features

* Monte Carlo–based QMRA framework
* Pathogen-specific dose–response modeling
* Configurable water sources, reuse activities, and treatment trains
* Explicit handling of degraded and failed treatment performance
* DALY-based health risk metrics (WHO-aligned)
* Interactive Streamlit interface
* CSV export and publication-ready figures
* Conservative aggregation of health risk across pathogens

---

## 🧪 Methodological Overview

The model estimates annual health risk per exposed individual by simulating:

1. Pathogen concentrations in source water
2. Exposure via selected reuse activities
3. Treatment performance (including degradation and failure states)
4. Infection and illness probabilities
5. Health burden expressed as **DALY per person-year**

### DALY Aggregation Across Pathogens

To avoid overestimation of health burden due to implicit assumptions of simultaneous disease episodes, **total DALY is aggregated using a max-per-iteration approach**:

> For each Monte Carlo iteration, the total DALY is defined as the maximum DALY contribution across all simulated pathogens.

This represents the most severe health outcome experienced by an individual in a given year and avoids unrealistic additive disease burden.

The summation of DALY across pathogens can be interpreted only as an upper-bound sensitivity case and is therefore not used as the default metric.

---

## 📊 Outputs

For each scenario, the model produces:

* Distribution of total DALY per person-year
* Median and 95th percentile DALY
* Comparison against WHO benchmark (1×10⁻⁶ DALY/person-year)
* Pathogen-wise DALY contributions
* Expected cases, infection probability, and annual infection risk
* Publication-ready plots
* Full Monte Carlo output as CSV

---

## 🖥️ User Interface

The Streamlit interface allows users to:

* Select water source(s)
* Select one or more reuse activities
* Select pathogens of concern
* Configure treatment trains and operational status
* Specify exposed population size
* Run simulations interactively
* Download results and figures

---

## 🔧 Installation

### Requirements

* Python ≥ 3.9
* Recommended: virtual environment

### Install dependencies

```bash
pip install -r requirements.txt
```

Key dependencies include:

* `numpy`
* `pandas`
* `scipy`
* `matplotlib`
* `streamlit`

---

## ▶️ Running the App

```bash
streamlit run app.py
```

The application will open in your browser.

---

## 📁 Repository Structure

```
QMRA_water_reuse/
│
├── app.py                             # Streamlit application
├── treatment_Scenarios_QMRA_*.py      # Core QMRA logic and parameters
├── requirements.txt
├── README.md
└── outputs/                           # (optional) exported results
```

---

## ⚠️ Intended Use and Limitations

* This tool is intended for **scenario analysis and comparative assessment**, not for regulatory approval without expert review.
* All risk calculations are deterministic given input assumptions; uncertainty is represented through Monte Carlo sampling.
* Results are sensitive to pathogen parameters, exposure assumptions, and treatment effectiveness distributions.

---

## 📚 References and Frameworks

* World Health Organization (WHO): *Guidelines for the Safe Use of Wastewater, Excreta and Greywater*
* Haas et al. (2014): *Quantitative Microbial Risk Assessment*
* ISO 16075 (Water reuse in irrigation)
* QMRA best practices in water reuse and drinking water safety planning

---

## 🤝 Contributing

Contributions are welcome, particularly in the form of:

* Additional pathogens or exposure pathways
* Alternative DALY aggregation approaches (clearly documented)
* Sensitivity and uncertainty analysis modules
* Validation against empirical data

Please open an issue or pull request to discuss proposed changes.

---

## 📄 License

MIT License

Copyright (c) 2026 Jesper Knutsson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


---

## 👤 Author

**Jesper Knutsson**
Researcher and consultant in sustainable water systems and water reuse

---

## 📝 Citation

If you use this tool in academic work, please cite as:

> Knutsson, J. (year). *QMRA Water Reuse Simulator*. GitHub repository: [https://github.com/jespknut/QMRA_water_reuse](https://github.com/jespknut/QMRA_water_reuse)

