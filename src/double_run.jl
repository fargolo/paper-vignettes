# Simultaneous test for recurrence and Poison times simultanously.
### To do: count number of recurrences that are part of diagonals or vertical sections

# To do: 1.  Iterate over diagonals 
#        2.  Fix indices for collumns (only upper half values)
#        3.  LinearAlgebra.diag(Matrix,n_diagonal)
"""
    double_inference(rec_matrix::RecurrenceMatrix;seqs="double")

Probabilities based on double runs.

Iterates over each column of a recurrence matrix and returns the probability of observing 
a consecutive run equal to or larger than the one in the data.
One can consider recurrences (seqs = "recurrences"), poincare times (seqs="poincare") 
or the largest among them (default: seqs = "double").  
"""

function double_inference(rec_matrix::RecurrenceMatrix;seqs="double")

    if seqs ∉ ["double","recurrences","poincare"]
        println("seqs must be either 'double', 'recurrences' or 'poincare'")
        return(NaN)
    end 
        
    mat_len = first(size(rec_matrix))
    probs = Float64[]
    

    # To do: implement map instead of for loop 
    # sequences = map(x-> diag(Matrix(rec_matrix),x), 1:size(rec_matrix)[1])

    for i in 1:(mat_len)
        cur_len = mat_len - i
        col_counts = StatsBase.rle(LinearAlgebra.diag(Matrix(rec_matrix),i))
        zipped_tups = collect(zip(col_counts[2],col_counts[1]))

        if seqs == "recurrences"
            zipped_tups = filter(i -> i[2] == 1, zipped_tups) # RECURRENCES: Tuples in which the 1st value is 1 
        elseif seqs == "poincare"
            zipped_tups = filter(i -> i[2] == 0, zipped_tups) # POINCARE TIMES: Tuples in which the 1st value is 1 
        end
        
        try
            max_val = maximum(first.(zipped_tups))
            cur_p = AnalyticComb.p_binary_words_doub_runl(max_val,cur_len)
            push!(probs,cur_p)
        catch e
            println("Sequence has zero occurences.")
            push!(probs,NaN)
        end         
        
    end

    return(probs)

end



function double_inference_weighted(rec_matrix::RecurrenceMatrix;seqs="double")

    if seqs ∉ ["double","recurrences","poincare"]
        println("seqs must be either 'double', 'recurrences' or 'poincare'")
        return(NaN)
    end 
        
    mat_len = first(size(rec_matrix))
    probs = Float64[]

    p = recurrencerate(rec_matrix)
    q = 1-p
    println("RR is ",p)
    println("P and Q are ",p," and ",q)

    

    # To do: implement map instead of for loop 
    # sequences = map(x-> diag(Matrix(rec_matrix),x), 1:size(rec_matrix)[1])

    for i in 1:(mat_len)
        cur_len = mat_len - i
        col_counts = StatsBase.rle(LinearAlgebra.diag(Matrix(rec_matrix),i))
        zipped_tups = collect(zip(col_counts[2],col_counts[1]))
        
        if seqs == "recurrences"
            zipped_tups = filter(i -> i[2] == 1, zipped_tups) # RECURRENCES: Tuples in which the 1st value is 1 
        elseif seqs == "poincare"
            zipped_tups = filter(i -> i[2] == 0, zipped_tups) # POINCARE TIMES: Tuples in which the 1st value is 1 
        end

        println("Zipped tuples: ")
        print(zipped_tups)
        println("Current segment length: ")
        print(cur_len)
        
        try
            max_val = maximum(first.(zipped_tups))
            cur_p = p_val_weighted(p,q,max_val,cur_len)
            push!(probs,cur_p)
        catch e
            println("Sequence has zero occurences.")
            push!(probs,NaN)
        end         
        
    end

    return(probs)

end




# Legacy
function double_inference_cols(rec_matrix::RecurrenceMatrix;seqs="double")

    if seqs ∉ ["double","recurrences","poincare"]
        println("seqs must be either 'double', 'recurrences' or 'poincare'")
        return(NaN)
    end 
        
    mat_len = first(size(rec_matrix))
    probs = Float64[]
    
    for i in 1:(mat_len)
        cur_len = mat_len - i + 1
        col_counts = StatsBase.rle(Matrix(rec_matrix)[1:cur_len,i])
        zipped_tups = collect(zip(col_counts[2],col_counts[1]))

        if seqs == "recurrences"
            zipped_tups = filter(i -> i[2] == 1, zipped_tups) # RECURRENCES: Tuples in which the 1st value is 1 
        elseif seqs == "poincare"
            zipped_tups = filter(i -> i[2] == 0, zipped_tups) # POINCARE TIMES: Tuples in which the 1st value is 1 
        end
        
        try
            max_val = maximum(first.(zipped_tups))
            cur_p = AnalyticComb.p_binary_words_doub_runl(max_val,cur_len)
            push!(probs,cur_p)
        catch e
            println("Sequence has zero occurences.")
            push!(probs,NaN)
        end         
        
    end

    return(probs)

end


function weighted_bin_runs_coeff(p,q,l,n)

    if (p+q - 1 > 0.01) || (l > n)
        println("p + q must be equal to 1 and l <= n")
            return(NaN)
    end 

    z = SymPy.symbols("z") 
    wei_mgf = (1 - p^l * z^l )/(1 - z + q*(p^l)*(z^(l+1))) # OGF
    coefs = collect(series(wei_mgf,z,0,n+1),z)
    coefs.coeff(z,n)

end


"""
    p_val_weighted(p,q,l,n)

p-value obtained from a one-tailed based on the exact distribution using the weighted model for consecutive runs `weighted_bin_runs_coeff`.  
"""
function p_val_weighted(p,q,l,n)

    if (p+q - 1 > 0.01) || (l > n)
        println("p + q must be equal to 1 and l <= n")
            return(NaN)
    end 

    weighted_coeffs = map(x->weighted_bin_runs_coeff(p,q,x,n),(l-1):1:n)
    probs = diff(weighted_coeffs)
    Float64(sum(probs))

end