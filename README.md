# A one-step spatial+ approach to mitigate spatial confounding in multivariate spatial areal models
This repository contains the R code to implement the methods described in the paper entitled "A one-step spatial+ approach to mitigate spatial confounding in multivariate spatial areal models" (Urdangarin et al., 2023) as well as the R code to create the figures and tables presented in the paper.


## Table of contents

- [Data](#Data)
- [Simulated data](#SimulatedData)
- [R code](#R-code)
- [References](#References)


# Data
### Rapes and dowry deaths data in Uttar Pradesh in 2011 [(Vicente et al., 2020)](https://rss.onlinelibrary.wiley.com/doi/10.1111/rssa.12545)

The [**data_UttarPradesh_2011.Rdata**](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/Data/data_UttarPradesh_2011.Rdata) file contains the following objects:
  - **_data_**: contains the data set used. It is a dataframe with the following variables:
    - **_ID_area_**: numeric identifiers of districts
    - **_dist_**: names of the districts of Uttar Pradesh
    - **_year_**: year in which the data is gathered
    - **_pop_**: population of each district in 2011
    - **_obs_**: number of rapes and dowry deaths in each district in 2011
    - **_exp_**: number of expected cases of rapes and dowry deaths in each district in 2011
    - **_Crime_**: 1=rapes, 2=dowry deaths
    - **_X1_**: standardized sex ratio covariate (number of females per 1000 males)
    - **_X5_**: standardized murder (per 100000 people) covariate 
    - **_X6_**: standardized burglary (per 100000 people) covariate 
    - **_X3_**: standardized female literacy rate (%) covariate

  - **_carto_UP_**: cartography of the 70 districts of Uttar Pradesh



# Simulated data
[Simulated_data](https://github.com/spatialstatisticsupna/Multivariate_confounding/tree/main/Simulated_data) folder contains a total of 11 .Rdata files (one file for each scenario) used in Simulation Study 1 and Simulation Study 2. Each .Rdata file contains the same objects as [**data_UttarPradesh_2011.Rdata**](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/Data/data_UttarPradesh_2011.Rdata) (**_data_** and **_carto_UP_**) but in Simulation study 1 a simulated covariate $X=(X_2, X_3)'$ is added to **_data_** whereas for Simulation study 2 the simulated covariate $X_1^{*}$ and the spatial effects $\theta=(\theta_1, \theta_2)'$ are added. Moreover, each .Rdata contains the following objects as well:

- **_log.risk_**: a vector that contains the simulated log risks for both crimes
- **_log.risk.crime1_**: a vector that contains the simulated log risks for crime 1
- **_log.risk.crime2_**: a vector that contains the simulated log risks for crime 2
- **_simu.O_**: a list with 300 simulated counts data sets for both crimes
- **_simu.O.crime1_**: a list with 300 simulated counts data sets for crime 1
- **_simu.O.crime2_**: a list with 300 simulated counts data sets for crime 2


# R code

R code to implement the M-models with one-step spatial+ procedure is available in the folder [**R**]. The folder contains the code to to fit all the models and reproduce the tables and figures of the paper. 

- [**R/Real_data_analysis**](https://github.com/spatialstatisticsupna/Multivariate_confounding/tree/main/R/Real_data_analysis) folder contains the R code used in the real data analysis.
  - [Figure1.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Real_data_analysis/Figure1.R): R script to reproduce Figure 1 of the paper.
  - [functions](https://github.com/spatialstatisticsupna/Multivariate_confounding/tree/main/R/Real_data_analysis/functions): contains the functions to fit the M-models implemented using rgeneric function of INLA.
  - [run_MICAR.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Real_data_analysis/run_MICAR.R): R script to fit the M-Spatial and M-SpatPlus models with ICAR prior.
  - [run_MPCAR.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Real_data_analysis/run_MPCAR.R): R script to fit the M-Spatial and M-SpatPlus models with PCAR prior.
  - [run_MBYM2.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Real_data_analysis/run_MBYM2.R): R script to fit the M-Spatial and M-SpatPlus models with BYM2 prior. 
  - [Tables_2_3_4.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Real_data_analysis/Tables_2_3_4.R): R code to reproduce Table 2, 3 and 4 of the paper.
 
- [**R/Simulation_Study_1**](https://github.com/spatialstatisticsupna/Multivariate_confounding/tree/main/R/Simulation_study_1) folder contains the R code used in Simulation Study 1. Before running the M-models, the option _Scenario_ (scenario 1, 2, 3, 4, 5 or 6) must be defined.
  - [SimuStudy1_simulate_data.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/SimuStudy1_simulate_data.R): R code to simulated the 300 counts datasets.
  - [Figure2.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/Figure2.R): R code to reproduce Figure 2 of the paper.
  - [run_MICAR_Spat.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/run_MICAR_Spat.R): R script to fit the M-Spatial model with ICAR prior to the 300 simulated datasets.
  - [run_MICAR_SpatPlus.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/run_MICAR_SpatPlus.R): R script to fit the M-SpatPlus models with ICAR prior to the 300 simulated datasets.
  - [run_MPCAR_Spat.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/run_MPCAR_Spat.R): R script to fit the M-Spatial model with PCAR prior to the 300 simulated datasets.
  - [run_MPCAR_SpatPlus.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/run_MPCAR_SpatPlus.R): R script to fit the M-SpatPlus models with PCAR prior to the 300 simulated datasets.
  - [run_MBYM2_Spat.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/run_MBYM2_Spat.R): R script to fit the M-Spatial model with BYM2 prior to the 300 simulated datasets.
  - [run_MBYM2_SpatPlus.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/run_MBYM2_SpatPlus.R): R script to fit the M-SpatPlus models with BYM2 prior to the 300 simulated datasets.
 
- 

  - [SimuStudy1_Null.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy1_Null.R): R script to fit the null model to the 100 simulated datasets.
  - [SimuStudy1_Spatial.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy1_Spatial.R): R script to fit the spatial model to the 100 simulated datasets.
  - [SimuStudy1_RSR.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy1_RSR.R): R script to fit the RSR model to the 100 simulated datasets.
  - [SimuStudy1_SpatialPlus_eigenvectors.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy1_SpatialPlus_eigenvectors.R): R script to fit SpatPlus5, SpatPlus10, SpatPlus15 and SpatPlus20 models to the 100 simulated datasets.
  - [SimuStudy1_SpatPlusP1.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy1_SpatPlusP1.R): R script to fit SpatPlusP1 model to the 100 simulated datasets.
  - [SimuStudy1_SpatPlusTP1.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy1_SpatPlusTP1.R): R script to fit SpatPlusTP1 model to the 100 simulated datasets.
  - [SimuStudy1_SpatPlusP2.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy1_SpatPlusP2.R): R script to fit SpatPlusP2 model to the 100 simulated datasets.
  - [SimuStudy1_SpatPlusTP2.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/SimuStudy1_SpatPlusTP2.R): R script to fit SpatPlusTP2 model to the 100 simulated datasets.
  - [Tables_6_7_8.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/Tables_6_7_8.R): R code to reproduce Table 6, 7 and 8 of the paper for each scenario and subscenario.
  - [Tables_supplementary_A1_A2_A3_A4.R](https://github.com/spatialstatisticsupna/Simulation_confounding_article/blob/main/R/Simulation_study_1/Tables_supplementary_A1_A2_A3_A4.R): R code to reproduce Tables A1, A2 and A3 of the supplementary material for each scenario and subscenario.
  


Computations were run using R-4.0.4, INLA version 21.02.23, mgcv version 1.8-40.

# Acknowledgements
This work has been supported by Project PID2020-113125RB-I00/ MCIN/ AEI/ 10.13039/501100011033.

![image](https://github.com/spatialstatisticsupna/Comparing-R-INLA-and-NIMBLE/blob/main/micin-aei.jpg)
 
# References

	
