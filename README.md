# CasTuner
Modelling of repression and derepression dynamics of CasTuner systems


To execute this code, run the following R files with the order below.
The raw data are also provided in locations specified in each R file.



Step 1a.
Estimate the dynamic of upregulation of Cas-Repressors upon dTAG-13 withdrawal.
Run file: " Fitting tagBFP upregulation.R  "

Step 1b.
Estimate the dynamic of degradation of Cas-Repressors upon dTAG-13 addition
Run file: " Fitting tagBFP downregulation.R  "

Step 2. 
To know what is the dose-response relationship between repressors and target gene at the steady-state (assumed to be day 4 of titration: data corresponding to Supplementary Fig. 3e)
Run file:  " Fitting hill curves to dose responses.R " 

Step 3.
Estimate mCherry degradation rate from CasRx time-course
Run file:  " ODE_REV.R" lines : 89-110

Step 4a.
Simulate ODE model using computed parameters from previous steps, to predict the dynamic of mCherry repression by each repressor system. Calculate if there is a delay from theoretical and simulated data.
Run file:  " ODE_KD.R " 

Step 4b.
Simulate ODE model using computed parameters from previous steps, to predict the dynamic of mCherry derepression by each repressor system. Calculate if there is a delay from theoretical and simulated data.
Run file:  " ODE_REV.R" 

Thanks, 
G.
