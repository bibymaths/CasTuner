# CasTuner

Code associated with the manuscript 

## CasTuner: a degron and CRISPR/Cas-based toolkit for analog tuning of endogenous gene expression
#### Gemma Noviello, Rutger A. F. Gjaltema and Edda G. Schulz  

This repository contains code and data for modelling of repression and derepression dynamics of CasTuner systems


To execute this code, run the following R files with the order below.
The raw data are also provided in locations specified in each R file.


Step 1a.
Estimate the dynamic of upregulation of Cas-Repressors upon dTAG-13 withdrawal.
Run file: " Fitting tagBFP upregulation.R  "

Step 1b.
Estimate the dynamic of degradation of Cas-Repressors upon dTAG-13 addition
Run file: " Fitting tagBFP downregulation.R  "

Step 1c. 
Fit the parameters of the dose-response relationship between repressors and target gene at the steady-state (assumed to be day 4 of titration: data corresponding to Supplementary Fig. 3e)
Run file:  " Fitting hill curves to dose responses.R " 

Step 2.
Simulate ODE model using computed parameters from previous steps, to predict the dynamic of mCherry derepression by each repressor system. This script also estimates the Cherry degradation rate from CasRx time-course. The script then tests whether there is a derepression delay by comparing experimental and simulated data.
Run file:  " ODE_REV.R" 

Step 3.
Simulate ODE model using computed parameters from previous steps, to predict the dynamic of mCherry repression by each repressor system. The script then tests whether there is a repression delay by comparing experimental and simulated data.
Run file:  " ODE_KD.R " 

The code will load the raw data from subfolder "fcs_files" and produce the plots found in the subfolder "plots" and the estimated parameters found in subfolder "parameters". Each script should run in < 5min on a standard desktop computer.


R, RStudio versions on which the code has been tested:

R version 3.6.3 (2020-02-29)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

R version 4.2.1 (2022-06-23)
Platform: x86_64-apple-darwin17.0     
Running under: macOS 12.5    

### Required packages
deSolve_1.30                
nlstools_2.0-0                         
wesanderson_0.3.6           
ggsci_2.9                  
ggridges_0.5.3              
RColorBrewer_1.1-2          
ggcyto_1.14.1               
flowWorkspace_3.34.1        
ncdfFlow_2.32.0            
openCyto_1.24.0             
flowCore_1.52.1             
extrafont_0.17             
nlme_3.1-152                
minpack.lm_1.2-1                                                 
egg_0.4.5                   
gridExtra_2.3                                           
dplyr_1.0.8                                                   
tidyr_1.1.3                
ggplot2_3.3.3               
tidyverse_1.3.1            


