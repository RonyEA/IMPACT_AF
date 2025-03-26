# IMPACTaf Simulation Model

The **IMPACTaf** model is an R-based, microsimulation and transition modeling framework designed to estimate individual-level health trajectories, simulate disease burden, and evaluate alternative health scenarios. It supports flexible configuration, synthetic population generation, and multi-country calibration using publicly available inputs.

---

## Installation

### 1. **Clone or Download the Repository**

```bash
git clone https://github.com/RonyEA/IMPACTaf.git
cd IMPACTaf
```

Or download as a ZIP and extract locally.

### 2. **Set Your Working Directory in R**

```r
setwd("/path/to/IMPACTaf")  # Replace with your path
```

### 3. **Run the Global Setup Script**

This installs dependencies, loads the simulation package, and builds documentation:

```r
source("global.R")
```

---

## How to Run a Simulation

### 1. Load Required Package

```r
library(IMPACTaf)
```

### 2. Initialize a Simulation

```r
IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")
```

The YAML file controls parameters like country, cohort size, simulation length, and aggregation.

### 3. Generate a Lifecourse

```r
IMPACTncd$
  del_logs()$
  del_outputs()$
  run(mc = 1:5, multicore = FALSE)
```

Output is saved in `./outputs/lifecourse/` as `.csv.gz`.

---

## Running Scenarios

You can simulate hypothetical interventions using `run_scenario()`:

```r
tm <- TransitionModel$new()

tm$run_scenario(
  prob_changes = list(AIS = 0.9),  # Reduce AIS probability by 10%
  scenario_name = "reduced_AIS_10pct"
)
```

Scenario outputs are saved under:

- `./outputs/transitions/scenarios/`
- `./outputs/lifecourse/transitions/scenarios/`

---

## 📁 Project Structure

```
IMPACTaf/
│
├── auxil/                  # Auxiliary utilities
├── backup/                 # Backup scripts (ignored from Git)
├── disease_burden/         # Country-wise disease burden inputs
├── inputs/                 # Base simulation inputs, including YAML config
├── outputs/                # Lifecourse and transition outputs
├── IMPACTaf_model_pkg/     # R package with simulation classes
│   └── R/                  # Class definitions: Simulation, Design, TransitionModel, etc.
├── vignettes/              # RMarkdown guides (setup, scenario simulation, plotting)
├── Scripts/                # Local scripts (ignored from Git)
├── global.R                # Setup script for installing and loading the package
├── README.md               # You're here!
```

---

## 📚 Vignettes and Guides

| Title | Description |
|-------|-------------|
| [Install Guide](https://ronyea.github.io/IMPACT_AF/how_to_install_model.html) | How to set up and run the model |
| [Run Lifecourse](https://ronyea.github.io/IMPACT_AF/how_to_create_lifecourse.html) | How to simulate health state trajectories |
| [Scenario Simulation](https://ronyea.github.io/IMPACT_AF/transitionmodel-overview.html) | Adjust transition probabilities and run what-if analyses |

---

## 🧾 License



---

## 🙌 Credits

Developed by [YOUR NAME]
If you use this model in research or policy analysis, please cite appropriately.
