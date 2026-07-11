# Modeling within-lineage evolution as a consequence of changes in the adaptive landscape across short and long timescales

__Authors:__ Marion Thaureau<sup>1†</sup> and Kjetil Lysne Voje<sup>1</sup>

__Affiliation:__ <sup>1</sup>Evolution and Paleobiology, Natural History Museum, University of Oslo

__Contact:__ <sup>†</sup>marion.thaureau@nhm.uio.no

__Article:__ [Thaureau and Voje, 2026](https://doi.org/10.1093/evolut/qpag127)

__Journal:__ Evolution

__Year:__ 2026  

__Abstract:__ The adaptive landscape has been suggested as a potential conceptual bridge between phenotypic evolution on generational to macroevolutionary timescales. However, this potential remains largely untapped due to a limited understanding of how the adaptive landscape changes across time. We assessed the dynamics of the adaptive landscape across different time intervals by analyzing phenotypic time series using multivariate models of adaptation. First, we examined whether a decrease in river water flow affected the optimal body mass of a salmon population over a few decades. Second, we explored whether temperature changes affected the optimal size of coccoliths in a species of coccolithophore across about a hundred thousand years in the Cretaceous. Finally, we analyzed the extent to which temperature changes affected the optimal shell size in a lineage of planktic foraminifera over a few million years during the Miocene. We find evidence that evolution induced by changes in the adaptive landscape, driven by environmental change, occurred in the salmon and foraminifer datasets, but we do not detect a causal link between the investigated paleo‑proxies and size changes in the coccolithophore lineage. We argue that applying the same modeling framework to datasets covering different time spans not only offers a promising way to better connect evolutionary processes across timescales but also enables tests of adaptive hypotheses about within-lineage evolution.

__Keywords:__ timeseries, Ornstein-Uhlenbeck process, adaptation-inertia framework, macroevolution, microevolution.

__Info:__ This repository contains scripts and data used for analyses in the publication.




## Description of the files:

Folders (*Data*, *Script*, *Results*) contain a subfolder for each analyzed dataset: the salmon dataset, the coccolithophore dataset, and the foraminifera dataset.



### 1 - *Data*

This folder contains the phenotypic and environmental time series and the evoTs objects.

* *Salmon\_data*

  * **Atlantic\_salmon\_size.txt** contains all salmon catch recorded in river Eira from 1925 to 2016, the year and exact date of the catch, and the body mass of the salmon in grams along with an ID number for each catch. Note that some ID numbers and exact dates of catch are not available, especially for the oldest data.
  * **Discharge.txt** is the recorded water flow volume in the Eira river from 1931 to 2016. It contains both a yearly average and an average from June to September in m³s⁻¹.
  * **salmon\_evoTS\_obj.RData** contains the multivariate evoTs object made of the salmon body mass and river water flow data, their variance, sample size, and time for each datapoint.

* *Coccolith\_data*

  * **Biscutum\_TotalcoccolithLength.txt** contains the coccolith length of the Cretaceous coccolithophore lineage *B. constans* measured in micrometers across 160,000 years. For every 50 datapoints in time, the file also contains the variance, the sample size (60), and the age of the sample.
  * **isotopes.txt** contains the geochemistry data measured for the same samples from which the coccolith length was measured. The file contains oxygen and carbon isotopic values in addition to CaCO3 and TOC values.
  * **coccolith\_evoTS\_obj.RData** contains the multivariate evoTs object made of the coccolith length and isotopic data, their variance, sample size, and time for each datapoint.

* *Foraminifera\_data*

  * **Hodell\_and\_Vayavananda\_1993\_Fohsella\_length.txt** contains the test length of the Miocene lineage *Globorotalia (Fohsella)* measured in micrometers across 3 million years. For every 114 datapoints in time, the file also contains the variance, the sample size (60), and the age of the sample.
  * **table\_d18O.csv** and **table\_d13C.csv** contain the isotopic values recorded in the same samples where the foraminifera test sizes were measured along with the variance, sample size (1; or 2 when averaged over two species), and age of the sample.
  * **foraminifera\_evoTS\_obj.RData** contains the multivariate evoTs object made of the foraminifera test length and isotopic data, their variance, sample size, and time for each datapoint.



### 2 - *Script*

This folder contains the R scripts used to run the different multivariate OU models. They all test either for independent evolution (both the trait and the environment are changing randomly); independent adaptation (the trait is evolving toward a fixed optimum); or dependent adaptation (the trait is evolving toward a moving optimum, and changes in that optimum are driven by the environmental variable.

The input variables for each model are the evoTs object containing the time series, and the matrices **A** and **Σ** (also called matrix **R** in the scripts), filled with either 0 or 1, which respectively define the absence/presence of causality and/or correlation.

Each model is run 100 times with different starting values to explore the loglikelihood landscape. The 100 iterations have sometimes been split into different scripts to improve efficiency while running the models. For example, in the scripts for model 2a, **salmon\_script\_model\_2a\_it1.R** contains iterations 1 to 50; and **salmon\_script\_model\_2a\_it50.R** contains iterations 50 to 100.



* *Salmon\_scripts*

  * **salmon\_script\_model\_\*.R** are the script files for each model run on the salmon dataset. We tested nine different models. One set of models allows both variables to follow an unbiased random walk (i.e., no deterministic part in the OU process). We ran versions without (model 1a) and with correlated random changes (model 1b) between the variables (i.e., non-zero off-diagonal elements in the Σ matrix (also called the R matrix in the scripts)). A second set of models allows the salmon body size to evolve as an OU process, while the water flow changes according to a random walk. We ran both versions with absence (model 2a) and presence of correlated random changes (model 2b), and a version where the water flow did affect the optimum of the salmon body size (dependent adaptation-model 4a).  Finally, a third set of models allows both salmon body size and water flow to vary as OU processes. Independent adaptation versions of this last set were run with the absence (model 3a) and presence of correlated changes (model 3b), as well as dependent adaptation versions where water flow affects the optimum of the salmon body mass without (model 5a) and with correlated random changes (model 5b).
  * **salmon\_script\_modelsp\_\*.R** and **salmon\_script\_modelsp2\_\*.R** are versions run in a second time, for which we set up additional input variables. We define some specific starting values obtained from the best model to increase the efficiency of the searching algorithm:  \_modelsp\_ is used when we parametrized the full **A** matrix and the diagonal elements of **Σ** matrix; \_modelsp2\_ is used when we parametrized only the diagonal elements of both **A** matrix and **Σ** matrix.
  * **salmon\_script\_results\_table.R** and **salmon\_script\_results\_table\_SE.R** contain the scripts used to produce the result tables for parameter estimates of the tested model, and their standard errors.



* *Coccolith\_scripts* and *Foraminifera\_scripts*

  * **coccolith\_script\_model\_\*.R** and **foraminifera\_script\_model\_\*.R** are the script files for each model, respectively, run on the coccolithophore and foraminifera datasets. We tested the same twelve models for both datasets. In the first set of models, all three variables (trait, oxygen isotope, and carbon isotope) change according to an unbiased random walk without (model 1a) and with correlated random stochastic changes (model 1b). In the second set of models, the size trait evolves as an OU process towards a fixed optimum (independent adaptation), while the environmental proxies change according to an unbiased random walk. In this set, we investigate both the absence (model 2a) and the presence of correlated random changes among all three variables (model 2b). The trait and the two environmental variables follow OU processes in the third set of models (Figure 4C). Two of the models in the third set describe independent adaptation in the trait without (model 3a) and with correlated random changes among the three variables (model 3b). Six of the models in the third set describe dependent adaptation without and with correlated random changes among all variables, where the oxygen isotopes (models 4a and 4b), carbon isotopes (models 5a and 5b), or both variables influence the trait optimum (models 6a and 6b)
  * **coccolith\_script\_modelsp2\_\*.R** and **foraminifera\_script\_modelsp2\_\*.R** are versions run in a second time, for which we set up starting values obtained from the best model, to increase efficiency of the searching algorithm. In these versions, we parametrized only the diagonal elements of both **A** matrix and **Σ** matrix.
  * **coccolith\_script\_results\_table.R**, **coccolith\_script\_results\_table\_SE.R**, **foraminifera\_script\_results\_table.R**, and **foraminifera\_script\_results\_table\_SE.R** contain the scripts used to produce the result tables for parameter estimates of the tested model, and their standard errors.



### 3 - *Results*

This folder contains the results for each model and tables summarizing the best parameter estimates with their standard errors.

The output variables include the loglikelihood and AICc, which allow us to compare the goodness of the model fit. We also obtain different parameter estimates characterizing each time series included in the model: z is the ancestral value, and θ is the primary optimum (for two time series such as the salmon body mass and the river waterflow, we get two z and two θ). Finally, we get output for the matrices **A** and **Σ** (dimensions corresponding to the number of time series in the model), which characterize the strength of causality and the strength of correlation between each pair of time series included in the modeling framework. A **A** matrix filled with 0 indicates an absence of causality and no dependent adaptation. A **Σ** filled with 0 indicates an absence of correlation.



* *Salmon\_results*

  * **salmon\_run\_model\_\*.txt** are the output text files, printed while the models were running. The time of convergence for each iteration and output for the best iteration are printed. This file is available for each model.
  * **salmon\_result\_table.csv** and **salmon\_result\_table\_SE.csv** are tables respectively giving the best parameter estimates and their standard errors for each model.
  * *S\_Rdatafiles* > **salmon\_allruns\_model\_\*.RData** contains the result output for each iteration. This file is available for each model.
  * *S\_Rdatafiles* > **salmon\_bestrun\_model\_\*.RData** contains the result output for the best iteration. This file is available for each model.
  * *S\_Rdatafiles* > **salmon\_bestmodel.RData** contains a duplicate of the best iteration for the best model output.



* *Coccolith\_results* and *Foraminifera\_results*

Results for the coccolithophore and foraminifera datasets are presented in the same way as results for the salmon dataset. The RData files are respectively called *C\_Rdatafiles* and *F\_Rdatafiles*.



### 4 - *Supplementary\_Material*

* *Simulation\_temporal\_trend*

  * **Simulation\_temporal\_trend.R** is the script simulating two independent time series sharing a temporal trend and then fitting different models (independent evolution, independent adaptation, dependent adaptation) to the simulated data (100 iterations).
  * **Results\_Simulation\_temporal\_trend.RData** is the result of the simulation saved in an R object.
  * **Results\_Simulation\_temporal\_trend.out** is the printed output of the results.



This folder includes a simulation demonstrating that the multivariate OU framework does not mistake a shared temporal trend for evidence of a causal relationship between variables.



* *Sensitivity\_analyses*

  * *Salmon\_SA*

    * *Data* > **salmon\_evoTS\_obj.RData** is the evoTs object built with the salmon body mass and Eira river water flow time series (value, variance, sample size, year).
    * *Data* > **salmon\_bestmodel.RData** is the best iteration for the best model output for the salmon dataset.
    * *Scripts* > **Simulation\_sensitivity\_salmon\_it\*R** is the script simulating time series with the output of the best model for the salmon dataset. It then fits models of independent evolution, independent adaptation, and dependent adaptation to the simulated data.
    * *Results* > **salmon\_sensitivity\_analysis\_it\*.RData** is the result for each iteration, saved as an R object. The 100 iterations are split into 4 different scripts.
    * *Results* > **result\_salmon\_sensitivity\_it\*.RData** is the printed output for each run.
  * *Foraminifera\_SA*

    * Results for the sensitivity analyses for the timeseries foraminifera test length - oxygen isotopic ratio are presented in the same way as results for the sensitivity analysis on the salmon dataset.



This folder contains the results for model recovery analyses we ran for the two cases where dependent adaptation was supported (salmon body mass – river water flow; foraminifera body length – oxygen isotope).





## Sharing/Access information

The data, scripts, and results are also available in a Dryad repository:

* \[https://doi.org/10.5061/dryad.95x69p91r]

For more information on how evoTS works:

* \[https://klvoje.github.io/evoTS/index.html]

Time series data were obtained from published studies:

* Salmon dataset: Jensen et al. (2022), \[https://doi.org/10.1073/pnas.2207634119]
* Coccolith dataset: Bornemann and Mutterlose (2006), \[https://doi.org/10.1016/j.geobios.2005.05.005]
* Isotopic data for the coccolith dataset: Bornemann et al. (2005), \[https://doi.org/10.1144/0016764903-171]
* Foraminifera dataset: Hodell and Vayavananda (1993), \[https://doi.org/10.1016/03778398(93)90019-T]

The time series can also be downloaded directly as an evoTs object from the Phenotypic Evolution Time Series database:

* \[https://pets.nhm.uio.no/PETS/]



## Code/Software

All scripts were run using R version 4.2.1 and the package evoTS version 1.0.3. Models were mostly run in parallel loops to optimize the efficiency of the algorithm. The computations were performed on resources provided by Sigma2 - the National Infrastructure for High-Performance Computing and Data Storage in Norway.
