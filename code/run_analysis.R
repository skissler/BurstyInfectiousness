library(tidyverse)										
library(odin)											

# Setup 
source('code/utils.R') 									
source('code/global_parameters.R')						
source('code/parameters.R') 							

# Psi inference
source('code/psi_inference.R') 							
source('code/psi_identifiability.R') 					
source('code/psi_identifiability_figures.R') 			

# Uncontrolled epidemics 
source('code/episims.R') 								
source('code/overdispersion_extinction.R') 				
source('code/overdispersion_heatmaps.R')				
source('code/survival.R')								

# Controlled epidemics 
source('code/isolation_te.R')							
source('code/isolation_overdispersion.R')				
source('code/isolation_gi_truncation.R')				
source('code/isolation_gi_truncation_deterministic.R')	
source('code/isolation_growth_rate.R') 					
source('code/gatheringsize_main.R') 					
source('code/gatheringsize_od.R') 						

# General parameter inference
source('code/growthrate.R') 							
source('code/g_identifiability.R')						

# Sensitivity analyses (robustness to burst model choice)
source('code/sensitivity_episims.R')
source('code/sensitivity_overdispersion.R')
source('code/sensitivity_isolation_te.R')
source('code/sensitivity_growthrate.R')
source('code/sensitivity_marginal_gi.R')
