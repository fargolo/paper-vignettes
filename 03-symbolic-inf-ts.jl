using Random, Distributions, RecurrenceAnalysis, SymbolicInference

#Hyperparams
rec_rate, sd_norm, p_window_size = 0.6, 0.3, 30
cycles = 12; 

# Time-series
linear_basis = (0:0.4:2*pi*cycles)
time_series_sin = map(sin, linear_basis)
n_sample = length(linear_basis)
dist = Distributions.Normal(0, sd_norm)
white_n_time_series = rand(dist, n_sample)
time_series_sin_noise = white_n_time_series .+ time_series_sin

# Rec mat
rec_mat = RecurrenceAnalysis.RecurrenceMatrix(time_series_sin_noise, 
                                rec_rate; fixedrate=true)

# Motifs inference
motifs_dict = SymbolicInference.rec_matrix_motifs(rec_mat; 
        window_range=collect(1:p_window_size), seqs="recurrences", n_motifs=2)
# Coordinates of motifs
coordinates_dict = extract_recurrences(time_series_sin_noise, motifs_dict; num_windows=30)
coordinates_dict[1]
p = SymbolicInference.plot_motifs(time_series_sin_noise, coordinates_dict; n_motifs=1)

# Persistence barcode
pl_persist_bar = persistence_barcode(time_series_sin_noise; n_windows=30)