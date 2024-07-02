using CSV 
using DataFrames
using SymbolicInference
using RecurrenceAnalysis
using Normalization
include("data_load.jl")

rec_rate = 0.3

normal_rec_matrix = RecurrenceAnalysis.RecurrenceMatrix(ts_raw_d, rec_rate; fixedrate=true)
ma_rec_matrix = RecurrenceAnalysis.RecurrenceMatrix(ts_ma_d, rec_rate; fixedrate=true)

# short range 
motifs_dict_normal = SymbolicInference.rec_matrix_motifs(normal_rec_matrix; window_range=collect(1:3), n_motifs=3)
motifs_dict_ma = SymbolicInference.rec_matrix_motifs(ma_rec_matrix; window_range=collect(1:3), n_motifs=3)


# long range 
motifs_dict_normal = SymbolicInference.rec_matrix_motifs(normal_rec_matrix; window_range=collect(1:12), n_motifs=2)
motifs_dict_ma = SymbolicInference.rec_matrix_motifs(ma_rec_matrix; window_range=collect(1:12), n_motifs=2)

coordinate_normal = SymbolicInference.extract_recurrences(ts_raw_d, motifs_dict_normal; num_windows=12)
coordinate_ma = SymbolicInference.extract_recurrences(ts_ma_d, motifs_dict_ma; num_windows=12)

p = SymbolicInference.plot_motifs(ts_raw_d, coordinate_normal; n_motifs=24)
p2 = SymbolicInference.plot_motifs(ts_ma_d, coordinate_ma; n_motifs=24)
