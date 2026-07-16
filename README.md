# How bursty infectiousness shapes epidemic dynamics
**[Stephen M. Kissler](mailto:stephen.kissler@colorado.edu)**

---

This repository contains the data and code for reproducing the findings from *How bursty infectiousness shapes epidemic dynamics and control*, by Stephen M. Kissler. 

## Package requirements
Code is fully implemented in R. The codebase requires the following packages: 

- **tidyverse** 
- **odin** 
- **igraph** 
- **parallel** 
- **patchwork** 
- **Rcpp** 
- **rstan** 
- **knitr** 

`Rcpp` is used to compile the C++ simulation kernel (`code/src/sim_stochastic_rcpp.cpp`, compiled automatically by `utils.R`), and `rstan` is used to fit the SARS-CoV-2 viral kinetics model (`code/sc2_kinetics_fit.stan`).

## Summary
The full set of analysis code is run by sourcing `code/run_analysis.R`. All code assumes it is being run from the project root (`BurstyInfectiousness/`). The code consists of the following sections: 

### Setup
Sources the project-specific utility functions and defines the global parameters.

- `utils.R` defines the project functions, and compiles the Rcpp simulation kernel (`src/sim_stochastic_rcpp.cpp`)
- `global_parameters.R` defines key simulation parameters (*e.g.,* default population size, number of simulations)
- `parameters.R` generates a list of attributes for the benchmark pathogens (influenza, SARS-CoV-2 omicron, measles) 

### Psi inference
Fits the SARS-CoV-2 viral kinetics model, then conducts the inference and identifiability analysis for the burstiness parameter $\psi$.

- `sc2_kinetics_fit.R` fits the SARS-CoV-2 viral kinetics model to the Ct-value data in `data/ct_dat_refined.csv` via a Stan MCMC (using `sc2_kinetics_fit.stan` and the helpers in `sc2_kinetics_utils.R`), and writes the posterior to `output/posterior_features_single_category.rds`. **This step is slow (tens of minutes to hours, several GB of RAM). It writes a large artifact that is not tracked in Git; once the `.rds` exists, the fit line can be commented out in `run_analysis.R` to skip refitting on subsequent runs.**
- `sc2_kinetics_figures.R` generates the viral kinetics figures from the saved posterior
- `psi_inference.R` conducts the inference of the burstiness parameter $\psi$ across ten pathogens using serial interval data from transmission clusters in the [OutbreakTrees](https://github.com/DrakeLab/taube-transmission-trees) database
- `psi_cluster_figures.R` assembles Figure 2 (panels A–E): the $\psi$ cluster inferences together with the viral kinetics infectiousness window
- `psi_identifiability.R` conducts analyses to determine the identifiability of $\psi$ from simulated data based on the OutbreakTrees dataset, including sensitivity to imperfect ascertainment and incubation period distribution misspecification
- `psi_identifiability_figures.R` generates figures for the $\psi$ identifiability analysis

The $\psi$ likelihood used by both `psi_inference.R` and `psi_identifiability.R` is defined in the shared helper `psi_likelihood.R`.

### Uncontrolled epidemics
Runs a basic set of uncontrolled epidemic simulations, calculates overdispersion with time-varying contacts, and conducts the survival analysis of time-to-epidemic-establishment.

- `episims.R` generates simulations of uncontrolled epidemics under the Gamma burst model for different values of $\psi$, and produces daily and cumulative case count figures
- `overdispersion_extinction.R` calculates the extinction probability of epidemics when contact rates vary over time, leading to coincidence superspreading
- `overdispersion_heatmaps.R` generates heatmaps of the simulation-based and analytic overdispersion resulting from time-varying contacts 
- `survival.R` conducts the epidemic establishment survival analysis, capturing the epidemic time shift that becomes more variable with bursty infectiousness. 

### Controlled epidemics
Runs a set of scripts to assess the impact of detect-and-isolate interventions and gathering size restrictions.

- `isolation_te.R` computes testing effectiveness for detect-and-isolate interventions with various detection strategies (deterministic, symptom-based, screening-based, and screening-based with turnaround delay) 
- `isolation_overdispersion.R` computes the overdispersion that results from symptom-based detection and isolation 
- `isolation_gi_truncation.R` computes the generation interval truncation that results from symptom-based detection and isolation using simulations of transmission clusters 
- `isolation_gi_truncation_deterministic.R` computes the generation interval truncation that results from deterministic detection and isolation analytically 
- `isolation_growth_rate.R` computes the adjusted growth rate that arises from generation interval truncation under D&I intervention
- `fig5_figures.R` assembles Figure 5 (panels B–E) from the saved detect-and-isolate panels
- `gatheringsize_main.R` generates simulations of gathering-size-restricted epidemics to compare extinction probabilities
- `gatheringsize_od.R` deterministically computes the overdispersion and predicted epidemic extinction probabilities produced by gathering size restrictions

### General parameter inference
Conducts simulation-based inference of the epidemic growth rate $r$ and the generation interval distribution $g(\tau)$.

- `growthrate.R` computes the epidemic growth rate $r$ under different levels of burstiness from epidemic simulations in an infinite population
- `g_identifiability.R` generates simulation-based inferences of the generation interval distribution parameters under diferent levels of burstiness

### Sensitivity analyses
Runs a set of sensitivity analyses to assess robustness to the choice of burst model.

- `sensitivity_episims.R` generates epidemic simulations under the alternative burst models
- `sensitivity_overdispersion.R` calculates overdispersion due to coincidence superspreading and produces overdispersion heatmaps under the alternative burst models 
- `sensitivity_isolation_te.R` calculates testing effectiveness under D&I interventions for the alternative burst models 
- `sensitivity_growthrate.R` estimates the epidemic growth rate $r$ under various levels of burstiness with the alternative burst models 
- `sensitivity_marginal_gi.R` compares transmission timing distributions for the alternative burst models to their marginal expectations to verify the stochastic simulation strategies 
- `ed_burstmodels_qq.R` builds the Extended Data QQ-plot grid from the marginal-GI objects produced above

### Extended Data assembly
Runs last, and must run after every analysis script above.

- `make_extended_data.R` takes the per-panel plot objects stashed by `save_fig()` as inputs, composites them into publication-formatted Extended Data figures under `figures/ExtendedData/`, and copies each into the manuscript figures directory

## A note on the figures
Sourcing `code/run_analysis.R` reproduces every analysis panel (via `save_fig()`) and all composited Extended Data figures (via `make_extended_data.R`). Figures 2 and 5 are composited by the called scripts `psi_cluster_figures.R` and `fig5_figures.R`. The final main-text Figures 1, 3, and 4 are hand-assembled outside R (Mathematica/Inkscape; e.g. `figures/gammapoisson.nb` generates Extended Data Fig. 3) and are not regenerated by sourcing `run_analysis.R`.

## Data
The `data/` directory contains the inputs required by the pipeline:

- `ct_dat_refined.csv` — SARS-CoV-2 Ct-value trajectories used to fit the viral kinetics model
- `outbreaktrees_data.RDS` / `data_tibble_trees.RDS` — transmission cluster data from the OutbreakTrees database used for the $\psi$ inference
