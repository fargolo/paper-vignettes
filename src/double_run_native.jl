

function double_inference_weighted(rec_matrix::RecurrenceMatrix;seqs="double",max_window=6)

    if seqs ∉ ["double","recurrences","poincare"]
        println("seqs must be either 'double', 'recurrences' or 'poincare'")
        return(NaN)
    end 
    

    mat_len = dim(rec_matrix)
    
    if mat_len < max_window
        println("max_window must be smaller than matrix length")
        return(NaN)
    end 
    

    
    probs = Float64[]

    p = recurrencerate(rec_matrix)
    q = 1-p
    println("RR is ",p)
    println("P and Q are ",p," and ",q)

    

    # To do: implement map instead of for loop 
    # sequences = map(x-> diag(Matrix(rec_matrix),x), 1:size(rec_matrix)[1])

    for i in 1:max_window
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
        println("\n Current segment length: ")
        print(cur_len)
        
        try
            max_val = maximum(first.(zipped_tups))
            println("Diagonal largest sequence size and total length are")
            print(max_val,cur_len)
            cur_p = p_val_weighted(p,q,max_val,cur_len)
            println("Probability:")
            print(cur_p)
            push!(probs,cur_p)
        catch e
            println("Sequence has zero occurences.")
            push!(probs,NaN)
        end         
        
    end

    return(probs)

end



function p_val_weighted(p,q,l,n)

    if (p+q - 1 > 0.01) || (l > n)
        println("p + q must be equal to 1 and l <= n")
            return(NaN)
    end 

    weighted_coeffs = map(x->weighted_bin_runs_coeff(p,q,x,n),(l-1):1:n)
    probs = diff(weighted_coeffs)
    Float64(sum(probs))

end

wei_mgf(p,q,l,z) = (1 - p^l * z^l )/(1 - z + q*(p^l)*(z^(l+1)))

function weighted_bin_runs_coeff(p,q,l,n)

    if (p+q - 1 > 0.01) || (l > n)
        println("p + q must be equal to 1 and l <= n")
            return(NaN)
    end 

    #z = TaylorSeries.set_variables("z") 
#   wei_mgf = (1 - p^l * z^l )/(1 - z + q*(p^l)*(z^(l+1))) # OGF
    #coefs = collect(series(wei_mgf,z,0,n+1),z)
    #getcoeff(taylor_expand(z -> wei_mgf(z),order=n_tot+1),n_tot)
    tay_exp = TaylorSeries.taylor_expand(z -> wei_mgf(p,q,l,z),order=n+1)
    TaylorSeries.getcoeff(tay_exp,n)

end

#@benchmark weighted_bin_runs_coeff(0.2,0.8,5,78)
#@benchmark weighted_bin_runs_coeff_bcp(0.2,0.8,5,78)

#function weighted_bin_runs_coeff_bcp(p,q,l,n)
#
#    if (p+q - 1 > 0.01) || (l > n)
#        println("p + q must be equal to 1 and l <= n")
#            return(NaN)
#    end 
#
#    z = SymPy.symbols("z") 
#    wei_mgf = (1 - p^l * z^l )/(1 - z + q*(p^l)*(z^(l+1))) # OGF
#    coefs = collect(series(wei_mgf,z,0,n+1),z)
#    coefs.coeff(z,n)
#
#end