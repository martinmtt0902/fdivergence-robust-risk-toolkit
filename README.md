Robust Risk Assessment Using F-Divergences

This repository contains a research-style project implementing a robust risk assessment toolkit based on F-divergences, inspired by:

Kruse, T., Schneider, L., & Schweizer, N. (2021).
“A Toolkit for Robust Risk Assessment Using F-Divergences.” Management Science.

The project includes theoretical study, simulation experiments, and an application to Hong Kong COVID-19 data.

Project Information

• Course: RMSC4102 – Research Project in Risk Management Science (CUHK)

• Type: Individual research project

• Language: R

• Focus: Robust risk assessment using F-divergence uncertainty sets


Project Overview

Project objectives:

• Study the F-divergence model uncertainty framework

• Implement key toolkit components (H, F, T functions, worst-case expectation, asymptotic variance)

• Reproduce and extend experiments from the paper

• Apply the approach to real data (COVID-19 HK), using a nonparametric bootstrap reference model

Main concepts:

• F-divergence ball around a reference model

• H-function inducing an F-divergence

• Worst-case expectation under model uncertainty

• Worst-case density described by α₁ and α₂


Repository Structure

Directory tree:
```text
├── code/
│  fdivergence_robust_risk_simulations.R
├── data/
│  COVIDHKdata.csv
├── report/
│  Project_Main.pdf
└── README.md
```

File descriptions:

• code/fdivergence_robust_risk_simulations.R — complete R implementation

• data/COVIDHKdata.csv — cleaned Hong Kong COVID-19 dataset

• report/Project_Main.pdf — full project report


Methods Implemented
1. Preliminary Functions

• alpha_est(): estimate worst-case density parameters

• WC_est(): compute worst-case expected value

• v_est(): asymptotic variance using the Delta method

• F_spline(): build F from H via numerical integration and spline interpolation

• Lognormal_Toolkit() and Weibull_Toolkit(): generate H, H⁻¹, F, T for lognormal and Weibull models


These form a compact F-divergence numerical toolkit.

2. Experiment 1 – Inventory Pooling (Lognormal)

• Lognormal demand model

• Compare decentralized vs centralized systems

• Compute nominal, worst-case, and best-case gains

• Add confidence bounds using asymptotic variance

• Plot of gains vs number of locations


3. Experiment 2 – Proportional Reinsurance under Dependence

• Dependence modeled with iid, Gaussian copula, and t-copula

• Weibull marginal distributions

• Compute worst-case expected reinsured losses for various κ

• Compare dependence structures


4. Additional Example 1 – Stability of Standardizers

• Poisson(π) sample

• Compare parametric vs non-parametric standardizers

• Compute worst-case expectations

• Sensitivity analysis under uncertainty


5. Additional Example 2 – Generalized Power Curve

• Mixture distribution as true model

• Compare KS test vs AD test

• Compute worst-case power under F-divergence balls

• Identify robustness differences


6. Additional Example 3 – COVID-19 Hong Kong Application

• Daily confirmed case counts (cleaned data)

• Nonparametric bootstrap reference model

• Compute worst-case expected mean and confidence band

• Real-data demonstration of the toolkit


Data Description

• COVIDHKdata.csv contains cleaned daily new confirmed cases in Hong Kong

• Single column: “x”

• Derived from Kaggle dataset “Novel Corona Virus 2019 Dataset” (by gpreda)


How to Run the Code

Install required packages:
```text
install.packages("nloptr")
install.packages("numDeriv")
install.packages("rpart")
```

Run the script:
```text
source("code/fdivergence_robust_risk_simulations.R")
```
This executes all experiments and outputs all figures.

Project Context and Contributions

• Implemented full F-divergence numerical toolkit

• Reproduced experiments from the reference paper

• Designed additional examples (standardizers, power curves, COVID HK application)

• Wrote the full report included in the repository


References

Kruse, T., Schneider, L., & Schweizer, N. (2021).
A Toolkit for Robust Risk Assessment Using F-Divergences. Management Science.

Kaggle: Novel Corona Virus 2019 Dataset (gpreda)
