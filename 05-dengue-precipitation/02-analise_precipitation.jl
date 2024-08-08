using CSV 
using DataFrames
using SymbolicInference
using RecurrenceAnalysis
using Normalization

rec_rate = 0.4;


SymbolicInference.persistence_motifs(ts_ma_max_p; range = collect(0.1:0.1:0.9), n_windows=52)
SymbolicInference.persistence_barcode(ts_ma_max_p; range = collect(0.1:0.1:0.9), n_windows=52,alpha_thresh=0.0001)


df = CSV.read("data/sanjuan_data.csv", DataFrame)
ts_raw = df[!, "precipitation"]
ts_ma = df[!, "mean_precip"]

ts_raw = convert(Vector{Float64}, ts_raw)
ts_ma = convert(Vector{Float64}, ts_ma[3:end])

norm_fit = fit(UnitEnergy, ts_raw)
ts_raw_unit = normalize(ts_raw, norm_fit)

norm_fit_ma = fit(UnitEnergy, ts_ma)
ts_ma_unit = normalize(ts_ma, norm_fit_ma)

# Matrix

normal_rec_matrix = RecurrenceAnalysis.RecurrenceMatrix(ts_raw, rec_rate; fixedrate=true)
ma_rec_matrix = RecurrenceAnalysis.RecurrenceMatrix(ts_ma, rec_rate; fixedrate=true)


motifs_dict_normal = SymbolicInference.rec_matrix_motifs(normal_rec_matrix; max_window=120, n_motifs=3)
motifs_dict_ma = SymbolicInference.rec_matrix_motifs(ma_rec_matrix; max_window=120, n_motifs=3)

coordinate_normal = SymbolicInference.extract_recurrences(ts_raw, motifs_dict_normal; num_windows=120)
coordinate_ma = SymbolicInference.extract_recurrences(ts_ma, motifs_dict_ma; num_windows=120)
# p = SymbolicInference.plot_motifs(ts_raw, coordinate_normal; n_motifs=4)
p2 = SymbolicInference.plot_motifs(ts_ma, coordinate_ma; n_motifs=6)
