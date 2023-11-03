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
