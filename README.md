# AgeStructuredRSVModelWithUrbanization
Code for "Urban contact patterns shape respiratory syncytial virus epidemics with implications for vaccination" by Kimball et al. Science Advances (2025, in press).


### There are five main rmd files to run here.

#### Data Analysis
This runs the basic data analysis done on the data and produces figures 1,2

#### SIR- Plots and example model run
Here, we run the SIR model and produce basic model output

#### SIR- Model Fit
This file is performs the model fitting using LHC random sampling to best fit Texas data

#### SIR- Bifurcation analysis
This file perturbs parameters and measure the effect on model output including periodicity and intensity

#### SIR- Vaccination analysis
Models the effect of vaccination on disease dynamics

### There are also 3 folders
#### output
stores data files that are outputed including bifurcation simulation data and fitted parameter data

#### source
contains scripts for all the functions referenced in the rmd files. See scripts for description of each.

#### data 
contains data files to run data analysis and SIR simulations

#### supplement-age_dist_figure
contains code to build age distribution of United States and ANOVA testing
