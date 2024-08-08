using TaylorSeries

function p_val_weighted(p,q,l,n)

    if (p+q - 1 > 0.01) || (l > n)
        println("p + q must be equal to 1 and l <= n")
            return(NaN)
    end 

    weighted_coeffs = map(x->weighted_bin_runs_coeff(p,q,x,n),(l-1):1:n)
    probs = diff(weighted_coeffs)
    Float64(sum(probs))

end

function weighted_bin_runs_coeff(p,q,l,n)

    if (p+q - 1 > 0.01) || (l > n)
        println("p + q must be equal to 1 and l <= n")
            return(NaN)
    end 

    wei_mgf(p,q,l,z) = (1 - p^l * z^l )/(1 - z + q*(p^l)*(z^(l+1)))
    tay_exp = TaylorSeries.taylor_expand(z -> wei_mgf(p,q,l,z),order=n+1)
    TaylorSeries.getcoeff(tay_exp,n)

end

function persistence_motifs(time_series; range = collect(0.1:0.1:0.9), n_windows=10)
    all_motifs = []

    for cur_range in range
        time_series_mat = RecurrenceAnalysis.RecurrenceMatrix(time_series, cur_range; fixedrate=true)
        motifs_result = SymbolicInference.rec_matrix_motifs(time_series_mat; max_window = n_windows, n_motifs=1)
        push!(all_motifs, motifs_result)
    end

    num_windows = length(all_motifs[1]["Probs"]) # Assuming all have the same number of motifs

    fig = Figure()
    ax = Axis(fig[1, 1], limits = (0, 1, 0, 1),
    xlabel="Recurrence Rates", 
        ylabel="Probabilities")

    for j in 1:num_windows
        x_vals = range
        y_vals = [all_motifs[i]["Probs"][j] for i in 1:length(range)]
        lines!(ax, x_vals, vcat(y_vals...), label="Window $j")
    end

    Legend(fig[1, 2], ax, "Motifs")
    fig
end
