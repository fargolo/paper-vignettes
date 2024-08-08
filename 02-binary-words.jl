using AnalyticComb
using Plots

# K occurrences
bin_words_with_k_occurences(4, 10)

# K occurrences contrs. d
bin_words_with_k_occurences_constr(4, 10, 2)

# Without k runs
words_without_k_run(9, 10; m=2)
counts = map(x -> 
    words_without_k_run(x, 10; m=2), 0:1:10)
probs = diff(counts) ./ 2^10
plot(probs)
sum(probs[7:10])


# Runs probs
probs = map(x -> 
    bin_words_runs_prob(x, 10), 1:1:10)

sum(probs[7:10])

# Weighted runs prob
raw_probs = map(x->
    weighted_bin_runs_prob(0.4, 0.6, x, 15), 0:1:15)
probs = diff(raw_probs)
plot(probs)
sum(probs[6:15])
sum(probs[7:15])
sum(probs[8:15])

weighted_bin_runs_pval(0.4, 0.6, 6, 15)
weighted_bin_runs_pval(0.4, 0.6, 7, 15)
weighted_bin_runs_pval(0.4, 0.6, 8, 15)
