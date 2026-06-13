# Code Repository for PhD Thesis

## Systemic Risk of Extreme Climate Events and Contagion of Economic Impacts

This repository contains all scripts used to produce the empirical analyses, simulations, figures, and results presented in the PhD thesis.

The objective of this research is to study how extreme climate events generate direct damages, propagate through economic networks, and contribute to systemic risk through global value chains.

---

# Repository Structure

## Chapter 2 – Direct Extreme Events Damage and Climate Models

### EM-DAT Statistical Analysis

Files in this section reproduce the descriptive analyses of historical disaster losses based on the EM-DAT database:

* descriptive statistics by disaster type
* temporal evolution of damages
* geographical distribution of losses
* damage distribution analysis
* preliminary risk indicators

### Extreme Value Theory

Scripts implementing:

* block maxima methods
* Generalized Extreme Value (GEV) estimation
* flood damage modelling
* tail risk estimation

### Climate Correlation Analysis

Scripts used to:

* compare damages with climate indicators
* analyze relationships with ERA5 variables
* project damage evolution under climate scenarios

---

## Chapter 3 – Indirect Extreme Events Damage

### Trade Network Construction

Scripts building international trade networks from trade-flow data.

### Network Metrics

Implementation of:

* centrality measures
* susceptibility indicators
* diffusion indicators
* network topology analysis

### Contagion Model

Simulation framework estimating indirect damages resulting from disruptions propagating through trade networks.

Outputs include:

* indirect losses
* contagion paths
* diffusion dynamics
* network vulnerability indicators

---

## Chapter 4 – Systemic Risk and Global Value Chains

### Global Value Chain Modelling

Scripts reproducing the conceptual and mathematical framework developed in Chapter 4.

### Network Simulations

Simulation experiments investigating:

* shock propagation
* amplification mechanisms
* buffering effects
* indirect loss generation

### Systemic Risk Indicators

Implementation of the proposed measures for assessing systemic exposure to climate-related disruptions.

---

## Chapter 5 – Climate Stress Testing

### Geographical Exposure Assessment

Scripts computing geographical exposure indicators using climate hazard information.

### Sectoral Exposure Assessment

Scripts estimating sector-specific vulnerability to climate extremes.

### Scenario Analysis

Implementation of stress scenarios based on:

* IPCC pathways
* NGFS scenarios
* CMIP climate projections

### Stress Testing Framework

Main implementation of the climate stress-testing model proposed in the thesis.

### ENGIE Case Study

Application of the framework to ENGIE assets, including:

* geographical exposure
* sectoral exposure
* stress calibration
* scenario results

---

# Reproducibility

The repository allows reproduction of:

* all descriptive statistics;
* all EVT analyses;
* all network indicators;
* all contagion simulations;
* all climate stress-test results;
* all figures and tables included in the thesis.

Input datasets were collected and pre-processed manually. Consequently, raw datasets are not necessarily included in this repository and may require separate access requests to their respective providers.

---

# Main Data Sources

* EM-DAT
* ERA5 Reanalysis
* CMIP Climate Projections
* BACI Trade Data
* UN Comtrade
* World Bank Climate Knowledge Portal
* NGFS Scenarios

---

# Author

Serine Guichoud

PhD Thesis

*Systemic Risk of Extreme Climate Events and Contagion of Economic Impacts*
