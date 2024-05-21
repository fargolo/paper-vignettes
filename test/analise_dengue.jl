using CSV 
using DataFrames
using SymbolicInference
using RecurrenceAnalysis
using Normalization

function extract_recurrences_v2(data_source::Vector{Float64}, 
    motifs_dict::Dict{String, Vector{Any}}; num_windows::Int64 = 3)
    
    probs_n = first(size(motifs_dict["Probs"][1]))
    probs_flat = vcat(motifs_dict["Probs"]...)
    sorted_probs_inds = sortperm(probs_flat)

    motifs_flat = vcat(motifs_dict["Motifs starts and duration"]...)
    motifs_by_window = motifs_flat[sorted_probs_inds[1:num_windows]]

    full_data_dict = []
    for (index, motif) in enumerate(motifs_by_window)
        motif_start, motif_size = motif
        motif_range = motif_start:(motif_start+motif_size)-1
        motif_y1 = data_source[motif_range]
        
        
        motif_window = div(sorted_probs_inds[index],probs_n)

        motif_rec_start = motif_start + motif_window
        motif_rec_range = motif_rec_start:(motif_rec_start+motif_size)-1 
        motif_y2 = data_source[motif_rec_range]
        
        push!(full_data_dict, Dict(
        "x1" => collect(motif_range), "y1" => motif_y1,
        "x2" => collect(motif_rec_range), "y2" => motif_y2,
        "window" => motif_window, "size" => motif_size,
        "prob" => probs_flat[sorted_probs_inds][index]))
    end
    return full_data_dict
end

rec_rate = 0.4;

df = CSV.read("data/sanjuan_data.csv", DataFrame)
ts_raw = df[!, "total_cases"]
ts_ma = df[!, "mean_cases"]

ts_raw = convert(Vector{Float64}, ts_raw)
ts_ma = convert(Vector{Float64}, ts_ma[3:end])

norm_fit = fit(UnitEnergy, ts_raw)
ts_raw_unit = normalize(ts_raw, norm_fit)

norm_fit_ma = fit(UnitEnergy, ts_ma)
ts_ma_unit = normalize(ts_ma, norm_fit_ma)

normal_rec_matrix = RecurrenceAnalysis.RecurrenceMatrix(ts_raw, rec_rate; fixedrate=true)
ma_rec_matrix = RecurrenceAnalysis.RecurrenceMatrix(ts_ma, rec_rate; fixedrate=true)

motifs_dict_normal = SymbolicInference.rec_matrix_motifs(normal_rec_matrix; max_window=120, n_motifs=3)
motifs_dict_ma = SymbolicInference.rec_matrix_motifs(ma_rec_matrix; max_window=120, n_motifs=3)

coordinate_normal = extract_recurrences_v2(ts_raw, motifs_dict_normal; num_windows=120)
coordinate_ma = extract_recurrences_v2(ts_ma, motifs_dict_ma; num_windows=120)
p = SymbolicInference.plot_motifs(ts_raw, coordinate_normal; n_motifs=4)
p2 = SymbolicInference.plot_motifs(ts_ma, coordinate_ma; n_motifs=6)
