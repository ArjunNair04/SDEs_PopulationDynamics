# Stochastic Differential Equations & Population Dynamics

Accompanying scripts for an undergraduate thesis that fit and diagnose stochastic differential‐equation (SDE) models to real population data. We cover:

1. **Microbial growth** under three media via a stochastic logistic SDE  
2. **Predator–prey dynamics** with the stochastic Lotka–Volterra (sLV) model  
3. **Predator–prey with saturation** via the stochastic Rosenzweig–MacArthur (sRM) model  

Each analysis script installs its own dependencies, reads the appropriate dataset, estimates drift & diffusion parameters, runs Euler–Maruyama simulations, and produces diagnostics.

---

## Repository Layout
SDEs_PopulationDynamics/
├── data/
│ ├── KHK growth curves_LB.xlsx
│ ├── KHK growth curves_MAA.xlsx
│ ├── KHK growth curves_M63.xlsx
│ └── Leigh1968_harelynx.csv
│
├── scripts/
│ ├── Bacteria_Modelling.R # stochastic logistic fit for E. coli OD data
│ ├── fit_lotka_volterra.R # stochastic Lotka–Volterra fit to lynx–hare
│ └── fit_rosenzweig_macarthur.R # stochastic Rosenzweig–MacArthur fit
│
├── output/
│ ├── plots/ # time series & simulation PNGs
│ └── tables/ # CSVs of parameter estimates & diagnostics
│
├── manuscript/ # LaTeX source of the write-up
└── README.md
