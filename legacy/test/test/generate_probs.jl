using Plots

# Time-series types:  
## Sine wave
## White noise
## Noisy sine wave
## Gaussian random walk (AR1)

sample_size = 300

series_probs = []

for my_seed in 1:sample_size
        # Hyperparams
        cycles = 8 # cycles on sin wave
        rec_rate = 0.2
        sd_norm = 0.5
        p_window_size = 50

        # Sine wave 
        time_series_sin = map(sin,0:0.5:2*pi*cycles)

        # White noise  
        rng = MersenneTwister(my_seed)
        n_sample=length(time_series_sin)
        dist = Normal(0,sd_norm)
        white_n_time_series = rand(dist,n_sample)

        # Noisy sine wave
        time_series_sin_nois = white_n_time_series .+ time_series_sin

        # Gaussian random walk (AR1)
        ## https://discourse.julialang.org/t/simple-random-walk-performance-tips/61553/2
        rand_walk_gauss_ts50 = rand_walk_gauss(n_sample;a=0.5,c=0.5)
        rand_walk_gauss_ts90 = rand_walk_gauss(n_sample;a=0.9,c=0.5)

        res_recs = map(x -> RecurrenceAnalysis.RecurrenceMatrix(x,rec_rate;fixedrate=true),
        [time_series_sin,white_n_time_series,time_series_sin_nois,rand_walk_gauss_ts50,rand_walk_gauss_ts90])
        all_probs = map(x-> double_inference_weighted(x;max_window=p_window_size,seqs="recurrences"),res_recs)
        push!(series_probs,all_probs)

end