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

### Uncontrolled epidemics
Runs a basic set of uncontrolled epidemic simulations, calculates overdispersion with time-varying contacts, and conducts the survival analysis of time-to-epidemic-establishment.

### Controlled epidemics
Runs a set of scripts to assess the impact of detect-and-isolate interventions and gathering size restrictions.

### Epi parameter inference
Conducts simulation-based inference of the epidemic growth rate $r$ and the generation interval distribution $g(\tau)$.

### Sensitivity analyses
*Runs a set of sensitivity analyses to assess robustiness to the burst model.* 