Code and data to replicate the results in "Robust Bayesian Nonparametric Variable Selection for Linear Regression".

This supplementary material consists of:

Source files
--------
BVVA.R (Spike-and-slab prior with Dirichlet mixtures)
hBVVA.R (Horseshoe prior with Dirichlet mixtures)
SVSS.R (Simple spike-and-slab)
hSVSS.R (Simple horseshoe)
BVVA_tstud.R (Spike-and-slab prior with Dirichlet mixtures with heavy tailed extension)
hBVVA_tstud.R (Horseshoe prior with Dirichlet mixtures with heavy tailed extension)
SVSS_tstud.R (Simple spike-and-slab with heavy tailed extension)
hSVSS_tstud.R (Simple horseshoe with heavy tailed extension)
GetSimulationData.R (Synthetic data generator)

Data files
--------
baseball.txt (Baseball dataset)
diabetes.tab.txt (Diabebtes dataset)
pollution.data.txt (Pollution dataset)
train.csv (Temperature dataset)
DREAM4 (Gene reconstruction dataset)
DREAM4gold (Gene reconstruction networks)

Implementation files
--------
synthetic_big.R (Scenario 1)
synthetic_big_var_3.R (Scenario 2)
k-foldCV.R (10-fold Cross Validation real world datasets)
k-foldCV_tstud.R (10-fold Cross Validation real world datasets with heavy tailed extension)
Source files to get every model's samples saved in folder "Results" along with the simulaton data.

DREAM4.R
Source file to get in working directory the samples in lists (network) of lists (gene) of lists named:
HSVS (Spike-and-slab prior with Dirichlet mixtures)
hHSVS (Horseshoe prior with Dirichlet mixtures)
SVS (Simple spike-and-slab)
hSVS (Simple horseshoe)
