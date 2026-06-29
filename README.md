# How bursty infectiousness shapes epidemic dynamics and control
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
- **knitr** 

## Summary
The full set of analysis code is run by sourcing `code/run_analysis.R`. All code assumes it is being run from the project root (`BurstyInfectiousness/`). The code consists of six main sections: 

### Setup
Sources the project-specific utility functions and defines the global parameters.

- `utils.R` defines the project functions 
- `global_parameters.R` defines key simulation parameters (*e.g.,* default population size, number of simulations)
- `parameters.R` generates a list of attributes for the benchmark pathogens (influenza, SARS-CoV-2 omicron, measles) 

### Psi inference
Conducts the inference and identifiability analysis for the burstiness parameter $\psi$.

- `psi_inference.R` conducts the inference of the burstiness parameter $\psi$ across ten pathogens using serial interval data from transmission clusters in the [OutbreakTrees](https://github.com/DrakeLab/taube-transmission-trees) database
- `psi_identifiability.R` conducts analyses to determine the identifiability of $\psi$ from simulated data based on the OutbreakTrees dataset, including sensitivity to imperfect ascertainment and incubation period distribution misspecification
- `psi_identifiability_figures.R` generates figures for the $\psi$ identifiability analysis

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
- `gatheringsize_main.R` generates simulations of gathering-size-restricted epidemics to compare extinction probabilities
- `gatheringsize_od.R` deterministically computes the overdispersion and predicted epidemic extinction probabilities produced by gathering size restrictions

### Epi parameter inference
Conducts simulation-based inference of the epidemic growth rate $r$ and the generation interval distribution $g(\tau)$.

- `growthrate.R` computes the epidemic growth rate $r$ under different levels of burstiness from epidemic simulations in an infinite population
- `g_identifiability.R` generates simulation-based inferences of the generation interval distribution parameters under diferent levels of burstiness

### Sensitivity analyses
*Runs a set of sensitivity analyses to assess robustiness to the burst model.* 

- `sensitivity_episims.R` Generates epidemic simulations under the alternative burst models
- `sensitivity_overdispersion.R` Calculates overdispersion due to coincidence superspreading and produces overdispersion heatmaps under the alternative burst models 
- `sensitivity_isolation_te.R` Calculates testing effectiveness under D&I interventions for the alternative burst models 
- `sensitivity_growthrate.R` Estimates the epidemic growth rate $r$ under various levels of burstiness with the alternative burst models. 
- `sensitivity_marginal_gi.R` Compares transmission timing distributions for the alternative burst models to their marginal expectations to verify the stochastic simulation strategies. 

