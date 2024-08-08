using Random, Distributions, RecurrenceAnalysis, SymbolicInference

rec_rate, sd_norm , p_window_size = 0.3 , 0.3 , 30
cycles = 12 ; 
linear_basis = (0:0.4:2*pi*cycles)
n_sample = length(linear_basis)
time_series_sin = map(sin,linear_basis)
rng = Random.MersenneTwister(42)
dist = Distributions.Normal(0,sd_norm)
white_n_time_series = rand(dist,n_sample)

time_series_sin_nois = white_n_time_series .+ time_series_sin

rec_mat = RecurrenceAnalysis.RecurrenceMatrix(
    time_series_sin_nois,rec_rate;fixedrate=true)

motifs_dict = SymbolicInference.rec_matrix_motifs(
    rec_mat;max_window=p_window_size,
    seqs="recurrences",n_motifs=2)

coordinates_dict = extract_recurrences(time_series_sin_nois, motifs_dict; num_windows=30)
p = SymbolicInference.plot_motifs(time_series_sin_nois, coordinates_dict; n_motifs=1)
using GLMakie
GLMakie.activate!()

p.center = false
GLMakie.save("test2.png",p)