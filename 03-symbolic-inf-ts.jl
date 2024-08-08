using Random, Distributions, RecurrenceAnalysis, SymbolicInference

rec_rate, sd_norm, p_window_size = 0.6, 0.3, 30
cycles = 12; linear_basis = (0:0.4:2*pi*cycles)
n_sample = length(linear_basis)
time_series_sin = map(sin, linear_basis)
rng = Random.MersenneTwister(42)
dist = Distributions.Normal(0, sd_norm)
white_n_time_series = rand(dist, n_sample)
time_series_sin_noise = white_n_time_series .+ time_series_sin