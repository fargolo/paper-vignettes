## A combinatorial multiverse of time-loops: Probabilistic inference on time-series based on recurrence analysis with symbolic methods

This repository contains code related to the [white-paper "A combinatorial multiverse of time-loops: Probabilistic inference on time-series based on recurrence analysis with symbolic methods"](https://osf.io/preprints/osf/3ws85).  

**Abstract**: Recurrence quantification analysis (RQA) is inspired by Poincaré's early studies and describes non-linear characteristics of dynamical systems from trajectories by identifying similarity between states pair-wisely among all observations. Since every pair is either close or distant, this procedure maps trajectories to the realm of binary states, a fruitful field for the application of combinatorial tools. 
In this work, we leverage symbolic methods from analytic combinatorics to make inference about the dynamics of systems represented as time-series data. Case studies include simulated data of Gaussian noise, periodic signal and autoregressive processes, as well as empirical data for precipitation volumes and dengue cases in San Juan, Puerto Rico. We demonstrate the detection of significant motifs: specific sequences of consecutive states that are repeated within a series or between two of them.
The framework successfully identifies patterns in systems such as random walks and noisy periodic signals. When applied to empirical data, it also highlights the association of dengue peaks (e.g., infectious outbreaks) and rainfall seasonality (e.g., weather fluctuation). 
Using combinatorial constructions tailored for special cases, our method provides exact probabilities for inference in the recurrence analysis of dynamical systems.  
Methods were implemented and made available in open-source software [AnalyticComb.jl](https://github.com/fargolo/AnalyticComb.jl) and [SymbolicInference.jl](https://github.com/fargolo/SymbolicInference.jl/).


