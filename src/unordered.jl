
function prob_occurrences_bin(k,n;alpha=0.05)
    probs = Float64[]
    counts = BigInt[]
    total_count = 2^n
    # Test: Sum probabilities of p(k), k > k_observed
    for i in 1:n
        cur_count = bin_words_with_k_occurences(i,n)
        cur_p = bin_words_with_k_occurences(i,n)/ total_count
#       print(Float64(cur_p)) 
        push!(probs,cur_p)
        push!(counts,cur_count)
    end
    return Dict("exact_prob" => sum(probs[k:n]), # exact probability
    "prob_distr" => probs, #prob. distribution
    "counts" => counts,
    "null_hypothesis" => ifelse(sum(probs[k:n]) < alpha,false,true)) #counts
end
