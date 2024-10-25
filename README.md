## A combinatorial multiverse of time-loops: Probabilistic inference on time-series based on recurrence analysis with symbolic methods

This repository contains code related to the [white-paper "A combinatorial multiverse of time-loops: Probabilistic inference on time-series based on recurrence analysis with symbolic methods"](https://osf.io/preprints/osf/3ws85).  

**Abstract**: Recurrence quantification analysis (RQA) is inspired by Poincaré's early studies and describes non-linear characteristics of dynamical systems from trajectories by identifying similarity between states pair-wisely among all observations. Since every pair is either close or distant, this procedure maps trajectories to the realm of binary states, a fruitful field for the application of combinatorial tools. 
In this work, we leverage symbolic methods from analytic combinatorics to make inference about the dynamics of systems represented as time-series data. Case studies include simulated data of Gaussian noise, periodic signal and autoregressive processes, as well as empirical data for precipitation volumes and dengue cases in San Juan, Puerto Rico. We demonstrate the detection of significant motifs: specific sequences of consecutive states that are repeated within a series or between two of them.
The framework successfully identifies patterns in systems such as random walks and noisy periodic signals. When applied to empirical data, it also highlights the association of dengue peaks (e.g., infectious outbreaks) and rainfall seasonality (e.g., weather fluctuation). 
Using combinatorial constructions tailored for special cases, our method provides exact probabilities for inference in the recurrence analysis of dynamical systems.  
Methods were implemented and made available in open-source software [AnalyticComb.jl](https://github.com/fargolo/AnalyticComb.jl) and [SymbolicInference.jl](https://github.com/fargolo/SymbolicInference.jl/).

# Notes from meeting - 19/09

-- Test method for logistic map, Roessler, Lorenz 63.  
    -- Compare with other methods (Symbolic dynamics?).
        -- E-mail: G. Corso (?)  

-- Paragraph on `run` vs. `sequence` (we use long-run to avoid ambiguity with sequence of SEQ operator).  

-- If mathematical technical specifications wont fit the text(e.g. notation, mid steps in algebraic formulas... ), add it to an Appendix.  

-- Provide more details about barcode plot
    -- Overall Rationale.  
    -- Colors in the plot (one for each window).   
    -- Extend text on how to interpret results of barcode plot for the example.  

-- Ask Roberto for commented version.  

-- See if there is a workaround for the duplicate figures including p-values. 
    -- There is no information about the thresholds in subtitle of the 1st figure. 
    -- It might be a good idea to include this information directly into the plot.   


-- When mentioning other operators, state that we will focus on the SEQ (Polya quasi-inverse) operator.
    -- There is broad field and the reader should refer to Flaj. and Sedg.'s book for more details and intuition about other operators (CYC, MSET, PSET)


-- Either insert reference to book equations or insert more steps about the derivation.   
 -- Break down algebraic expressions such as the one in Equation 3.
 -- State that atomic symbols will have arbritrary size "z"

-- Aim of paper: State in the beggining that we are connecting rec. matrices and analytic combinatorics.  
-- Extend explanations about the examples.    
-- Definition of a generating function.  
-- Suani topic 5) About the requirement  P.3, could you explain what the operator SEQ is and how SEQ is defined and used more deeply?It is important to present formally the  Definition of SEQ in the last paragraph of section III."  
-- -- The notion of GF is also mentioned in the Euler solution to polygon triangs.  
-- -- This section may also be a good one to introduce the concept of GF.  

-- Improve figure quality and interpretation of them.  
    -- Figure/label size.  

#### Notas - Reuniao

Incluir bordas (?)

Fig. 1 - Measure x Time - Explicar o que eh o eixo y
 --- Considerar explicar no label. The measure standards for simulated data that could apply to any time-series dataset.
-- Mudar tamanho da fonte

 
Fig. 2 -  Aumentar o tamanho da pagina 2
  -- Considerar remover.

Fig. 3 - 

Fig 4, 5, 6 - Eixo 

Fig. 9 - Adicionar Gaussiana 

Fig. 8, 10, 11 - Figuras em tamanhos diferentes

Windows - Definicao melhor de time-window (lag t-1, t-2, ..., t-k).

--

