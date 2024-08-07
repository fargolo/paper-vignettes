"""
persistence_barcode(time_series1,time_series2; range = collect(0.1:0.1:0.9), n_windows=10,alpha_thresh=0.05)

Return barcode plot for each window. 
    Points are plotted whenever the p-value is smaller than `alpha_thresh`.  

"""
function persistence_barcodev2(time_series1,time_series2; 
    range = collect(0.1:0.1:0.9), n_windows=10,alpha_thresh=0.05)
    all_motifs = []

    for cur_range in range
        time_series_mat = JointRecurrenceMatrix(time_series1,time_series2, cur_range; fixedrate=true)
        motifs_result = rec_matrix_motifs(time_series_mat; window_range=collect(1:n_windows), n_motifs=1)
        push!(all_motifs, motifs_result)
    end

    fig = Figure()
    ax = Axis(fig[1, 1], limits=(0, 1, 0, n_windows + 1),
              xlabel="Recurrence Rates", ylabel="Windows")

    for cur_window in 1:n_windows
        x_vals = range
        y_vals = [all_motifs[j]["Probs"][cur_window] for j in 1:length(range)]
        inds = findall(x -> x < alpha_thresh, vcat(y_vals...))
        scatter!(ax, x_vals[inds], fill(cur_window,length(inds)) , label="Window $cur_window")

    end
    fig
end
