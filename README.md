# SymbolicInference.jl
Probability-based inferences based on the symbolic method.  

White-paper link: https://docs.google.com/document/d/1giF2vcG8dcEHf4SONktg21Kk9KU5SXEAyO0fpNnyINY/edit#heading=h.d898sqp17elx  

## Basic pipeline

1 - Assume a probability distribution for recurrences and non-recurrences. Currently, only the uniform distribution ($RR = 0.5$, $P(1)=P(0)$ or $p = q$) is implemented.  
2 - Given the overall probabilities ($p , q$) of recurrences, check for underlying autocorrelation patterns by making inference over the average and maximum lenght of consecutive runs on diagonals.  
3 - Fun and profit.  

## To do (Roadmap)  

1 - Lower-order windows only?  

2 - Results on simulated data.  
-- Gaussian Random Walk  
-- White Noise
-- Periodic signal (sine wave)
-- Noisy periodic signal (white noise + sine wave)  

The `double_inference` function uses the solution for the consecutive runs problem (Flaj. & Sedgewick, page 52.   and R ÉV ÉSZ , P. Strong theorems on coin tossing. In Proceedings of the International Congress of Mathematicians (Helsinki, 1978) (Helsinki, 1980), Acad. Sci. Fennica, pp. 749–754.)  

Run the `test_paral_weighted.jl` file to obtain current results.  