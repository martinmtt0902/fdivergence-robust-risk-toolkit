# F-Divergence Robust Risk Assessment Toolkit

This repository contains an implementation of a **robust risk assessment framework** based on **F-divergences**, inspired by the paper  
**“A Toolkit for Robust Risk Assessment Using F-Divergences” (Kruse, Schneider & Schweizer, Management Science 2021).**

The project reproduces key numerical experiments from the paper and extends them with additional simulation studies and a real-data analysis using Hong Kong COVID-19 case counts.  
All computations are implemented in **R**.

The aim of this repository is to demonstrate practical understanding of divergence-robust optimization, stress testing under model uncertainty, and robust statistical inference using transform-based techniques.

---

## Project Information

- **Completion date:** May 2023  
- **Course:** RMSC4102 – Research Project
- **Languages:** R

---

## Overview

The project implements a modular toolkit for:

- Constructing **F-divergences** from an associated convex generator  
- Computing **worst-case expectations** under distributional uncertainty  
- Estimating **density tilting parameters** via constrained nonlinear optimization  
- Computing **asymptotic variances** using matrix Delta methods  
- Applying the framework to:
  - **Inventory pooling** under lognormal demand  
  - **Proportional reinsurance** with different copula dependencies  
  - **Robust comparison of test powers** (KS vs AD)  
  - **Empirical robust inference** for COVID-19 data

The implementation follows the structure of the original paper but is fully reproduced in R without external dependencies besides `nloptr` and `numDeriv`.

---

## Repository Structure

fdivergence-robust-risk-toolkit/
├── code/
│ └── fdivergence_robust_risk_simulations.R
│
├── data/
│ └── COVIDHKdata.csv
│
├── report/
│ └── Project_Main.pdf
│
└── README.md


### Explanation of Files

- **code/fdivergence_robust_risk_simulations.R**  
  Main R script implementing the full toolkit, simulations, and real-data application.

- **data/COVIDHKdata.csv**  
  Cleaned daily confirmed COVID-19 cases for Hong Kong (derived from Kaggle dataset).

- **report/Project_Main.pdf**  
  Full project report: theoretical summary, methodology, experiment replication, additional examples, and conclusions.

---

## Methods Implemented

### 1. Preliminary Functions

The script constructs a complete functional toolkit to perform F-divergence computations:

- **alpha_est()**  
  Estimates the worst-case tilting parameters  
  \(\alpha = (\alpha_1, \alpha_2)\)  
  by solving a constrained optimization problem using `nloptr`.

- **WC_est()**  
  Computes the **worst-case expectation** of a random variable under the tilted density.

- **v_est()**  
  Computes the **asymptotic variance** of the worst-case estimator using a matrix Delta method.

- **F_spline()**  
  Numerically constructs an **F-divergence** from a generator \( H \) via spline integration.

- **Lognormal_Toolkit()** and **Weibull_Toolkit()**  
  Provide the functions \(H\), \(H^{-1}\), \(F\), and the associated transform \(T\) for lognormal and Weibull models.

These functions collectively form a **compact and reusable F-divergence toolkit**.

---

## 2. Experiment 1 – Inventory Pooling (Lognormal Case)

This experiment reproduces the inventory pooling example from the paper.

- Lognormal demand model  
- Decentralized vs centralized inventory systems  
- Computes nominal, worst-case, and best-case pooling gains  
- Constructs confidence bounds using asymptotic variance  
- Reproduces Figures similar to those in Kruse et al. (2021)

---

## 3. Experiment 2 – Proportional Reinsurance Under Copula Dependence

Simulation of proportional reinsurance losses with:

- iid claims  
- Gaussian copula dependence  
- t-copula dependence  

For each dependence structure:

- Generate multivariate claim vectors  
- Standardize using mean/variance  
- Compute robust worst-case expectations over a grid of \(\kappa\) divergence radii  
- Compare sensitivity across dependence regimes

---

## 4. Additional Example 1 – Stability of Standardizers

A methodology diagnostic: assess robustness of two standardization methods

- Parametric: \(\sqrt{E[X]}\)  
- Non-parametric: \(\sqrt{\mathrm{Var}(X)}\)

The worst-case expectations reveal the sensitivity of each approach to heavy-tailed perturbations.

---

## 5. Additional Example 2 – Generalized Power Curve (Test Power Comparison)

Simulation comparing the **KS-test** and **AD-test** under mixture alternatives.

- Generates mixture normal data  
- Computes empirical rejection indicators  
- Applies F-divergence robustification to the empirical power  
- Produces worst-case power curves for each test  
- Shows which test is more robust to distributional misspecification

---

## 6. Additional Example 3 – COVID-19 Data Application

Real-data demonstration using Hong Kong daily confirmed COVID-19 cases.

- Bootstrap sampling to form empirical distribution  
- Compute robust worst-case mean estimates  
- Construct confidence intervals using asymptotic variance  
- Visualize sensitivity of mean case counts to distributional uncertainty  
- Shows practical value of F-divergence robustification in epidemiological risk assessment

---

## Purpose of This Repository

This repository demonstrates my ability to:

- Implement advanced risk analytics tools in R  
- Reproduce research-level numerical studies  
- Understand robust optimization and model uncertainty  
- Apply divergence-based techniques to both simulation and real-world data  
- Write clean, modular, and well-documented statistical code  

These skills directly translate to data science, quantitative research, actuarial analytics, and risk management roles.

---

## Disclaimer

The COVID-19 data file included here is a **cleaned subset** of the public Kaggle dataset:  
https://www.kaggle.com/datasets/gpreda/coronavirus-2019ncov  
As required by Kaggle’s license, the full dataset is **not** redistributed.

All code and materials are for academic and educational use only.

