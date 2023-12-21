# A simplified spatial+ approach to mitigate spatial confounding in multivariate spatial areal models
This repository contains the R code to implement the methods described in the paper entitled "A simplified spatial+ approach to mitigate spatial confounding in multivariate spatial areal models" (Urdangarin et al., 2023) as well as the R code to create the figures and tables presented in the paper.


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
The [**Simulated_data**](https://github.com/spatialstatisticsupna/Multivariate_confounding/tree/main/Simulated_data) encompasses a collection of 11 .Rdata files, with each file corresponding to a distinct scenario employed in Simulation Study 1 and Simulation Study 2. Each .Rdata file contains the same objects as [**data_UttarPradesh_2011.Rdata**](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/Data/data_UttarPradesh_2011.Rdata) (**_data_** and **_carto_UP_**). However, in Simulation study 1, a simulated covariate $X=(X_2, X_3)'$ is added to **_data_**, whereas Simulation Study 2 involves the inclusion of both the simulated covariate $X_1^{*}$ and the spatial effects $\theta=(\theta_1, \theta_2)'$. Additionally, every .Rdata file also accommodates the subsequent objects:
- **_log.risk_**: a vector that contains the simulated log risks for both crimes
- **_log.risk.crime1_**: a vector that contains the simulated log risks for crime 1
- **_log.risk.crime2_**: a vector that contains the simulated log risks for crime 2
- **_simu.O_**: a list with 300 simulated counts data sets for both crimes
- **_simu.O.crime1_**: a list with 300 simulated counts data sets for crime 1
- **_simu.O.crime2_**: a list with 300 simulated counts data sets for crime 2


# R code

The folder labeled [**R**](https://github.com/spatialstatisticsupna/Multivariate_confounding/tree/main/R) holds the necessary R code for executing the M-models using the simplified spatial+ approach. The folder contains the code to fit all the models and reproduce the tables and figures of the paper. 

- [**R/Real_data_analysis**](https://github.com/spatialstatisticsupna/Multivariate_confounding/tree/main/R/Real_data_analysis) folder contains the R code used in the real data analysis.
  
  The main files to fit the M-Spatial and M-SpatPlus models with ICAR, PCAR and BYM2 priors are [run_MICAR.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Real_data_analysis/run_MICAR.R), [run_MPCAR.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Real_data_analysis/run_MPCAR.R) and [run_MBYM2.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Real_data_analysis/run_MBYM2.R) respectively.
  - [Figure1.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Real_data_analysis/Figure1.R): R script to reproduce Figure 1 of the paper.
  - [functions](https://github.com/spatialstatisticsupna/Multivariate_confounding/tree/main/R/Real_data_analysis/functions): folder that contains the functions of M-models implemented using rgeneric function of INLA.
  - [Tables_2_3_4.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Real_data_analysis/Tables_2_3_4.R): R code to reproduce Table 2, 3 and 4 of the paper. Before running the models, the `spatial` argument (one of either "ICAR", "PCAR" or "BYM2") must be defined at the top of the code.
 
- [**R/Simulation_Study_1**](https://github.com/spatialstatisticsupna/Multivariate_confounding/tree/main/R/Simulation_study_1) folder comprises the R code employed during Simulation Study 1. Before running the scripts, the `scenario` argument (one of either "Scenario1", "Scenario2", "Scenario3", "Scenario4", "Scenario5" or "Scenario6") must be defined.
  - [SimuStudy1_simulate_data.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/SimuStudy1_simulate_data.R): R code to simulate the 300 counts datasets for crime 1 and crime 2.
  - [Figure2.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/Figure2.R): R code to reproduce Figure 2 of the paper.
  - [run_MICAR_Spat.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/run_MICAR_Spat.R): R script to fit the M-Spatial model with ICAR prior to the 300 simulated datasets.
  - [run_MICAR_SpatPlus.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/run_MICAR_SpatPlus.R): R script to fit the M-SpatPlus models with ICAR prior to the 300 simulated datasets.
  - [run_MPCAR_Spat.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/run_MPCAR_Spat.R): R script to fit the M-Spatial model with PCAR prior to the 300 simulated datasets.
  - [run_MPCAR_SpatPlus.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/run_MPCAR_SpatPlus.R): R script to fit the M-SpatPlus models with PCAR prior to the 300 simulated datasets.
  - [run_MBYM2_Spat.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/run_MBYM2_Spat.R): R script to fit the M-Spatial model with BYM2 prior to the 300 simulated datasets.
  - [run_MBYM2_SpatPlus.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/run_MBYM2_SpatPlus.R): R script to fit the M-SpatPlus models with BYM2 prior to the 300 simulated datasets.
  - [SimuStudy1_merge_results.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/SimuStudy1_merge_results.R): R code to combine the models fitted across 300 simulated datasets into a single list.
  - [Tables_5_6_7_8_9.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/Tables_5_6_7_8_9.R): R code to reproduce Table 5, 6, 7, 8 and 9 of the paper for each scenario and prior. Before running the code, the `spatial` argument (one of either "MICAR", "MPCAR" or "MBYM2") must be defined at the top of the code.
  - [Figure4.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/Figure4.R): R code to reproduce the boxplots in Figure 4.
  - [Figure5.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/Figure5.R): R code to reproduce the boxplots in Figure 5.
  - [FigureA1_supplementary.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/FigureA1_supplementary.R): R code to reproduce Figure A.1 in the supplementary material.
  - [Tables_A1toA11_supplementary.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/Tables_A1toA11_supplementary.R): R code to reproduce Tables A.1 to A.11 of the supplementary material A for each scenario and prior.
  - [FigureA2_supplementary.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/FigureA2_supplementary.R): R code to reproduce the boxplots in Figure A.2.
  - [FigureA3_supplementary.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_1/FigureA3_supplementary.R): R code to reproduce the boxplots in Figure A.3.
 
- [**R/Simulation_Study_2**](https://github.com/spatialstatisticsupna/Multivariate_confounding/tree/main/R/Simulation_study_2) folder comprises the R code employed during Simulation Study 2. Before running the scripts, the `scenario` argument (one of either "Scenario1", "Scenario2", "Scenario3", "Scenario4" or "Scenario5") must be defined.
  - [SimuStudy2_simulate_data.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_2/SimuStudy2_simulate_data.R): R code to simulate the 300 counts datasets for crime 1 and crime 2.
  - [Figure3.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_2/Figure3.R): R code to reproduce Figure 3 of the paper.
  - [run_MICAR_Spat.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_2/run_MICAR_Spat.R): R script to fit the M-Spatial model with ICAR prior to the 300 simulated datasets.
  - [run_MICAR_SpatPlus.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_2/run_MICAR_SpatPlus.R): R script to fit the M-SpatPlus models with ICAR prior to the 300 simulated datasets.
  - [run_MPCAR_Spat.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_2/run_MPCAR_Spat.R): R script to fit the M-Spatial model with PCAR prior to the 300 simulated datasets.
  - [run_MPCAR_SpatPlus.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_2/run_MPCAR_SpatPlus.R): R script to fit the M-SpatPlus models with PCAR prior to the 300 simulated datasets.
  - [run_MBYM2_Spat.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_2/run_MBYM2_Spat.R): R script to fit the M-Spatial model with BYM2 prior to the 300 simulated datasets.
  - [run_MBYM2_SpatPlus.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_2/run_MBYM2_SpatPlus.R): R script to fit the M-SpatPlus models with BYM2 prior to the 300 simulated datasets.
  - [SimuStudy2_merge_results.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_2/SimuStudy2_merge_results.R): R code to combine the models fitted across 300 simulated datasets into a single list.
  - [Tables_10_11_12_13_14.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_2/Tables_10_11_12_13_14.R): R code to reproduce Table 10, 11, 12, 13 and 14 of the paper for each scenario and prior. Before running the code, the `spatial` argument (one of either "MICAR", "MPCAR" or "MBYM2") must be defined at the top of the code.
  - [Figure6.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_2/Figure6.R): R code to reproduce the boxplots in Figure 6.
  - [Figure7.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_2/Figure7.R): R code to reproduce the boxplots in Figure 7.
  - [FigureB4_supplementary.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_2/FigureB4_supplementary.R): R code to reproduce Figure B.4 in the supplementary material.
  - [Tables_B12toB22_supplementary.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_2/Tables_B12toB22_supplementary.R): R code to reproduce Tables B.12 to B.22 of the supplementary material B for each scenario and prior. Before running the code, the `spatial` argument (one of either "MICAR", "MPCAR" or "MBYM2") must be defined at the top of the code.
  - [FigureB5_supplementary.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_2/FigureB5_supplementary.R): R code to reproduce the boxplots in Figure B.5.
  - [FigureB6_supplementary.R](https://github.com/spatialstatisticsupna/Multivariate_confounding/blob/main/R/Simulation_study_2/FigureB6_supplementary.R): R code to reproduce the boxplots in Figure B.6.
 


Computations were run using R-4.2.1, INLA version 22.12.16 (dated 2022-12-23).

# Acknowledgements
This work has been supported by Project PID2020-113125RB-I00/ MCIN/ AEI/ 10.13039/501100011033.

![image](https://github.com/spatialstatisticsupna/Comparing-R-INLA-and-NIMBLE/blob/main/micin-aei.jpg)
 
# References
[Urdangarin, A., Goicoa, T. , Kneib, T. and Ugarte, M.D. (2023). A simplified spatial+ approach to mitigate spatial confounding in multivariate spatial areal models. _Spatial Statistics (in press)_, DOI: 10.48550/arXiv.2308.11260.](https://doi.org/10.48550/arXiv.2308.11260)
	
