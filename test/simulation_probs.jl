using Plots

using Pkg
Pkg.add(path="/home/boitata/sandbox/SymbolicInference.jl")
using SymbolicInference

include("/home/boitata/sandbox/paper-vignettes/src/VignetteDependencies.jl")

# Script hyperparams
sample_sizes = [300,1000]
all_cycles = [4,8,16] # cycles on sin wave

for  cycles in all_cycles
        for sample_size in sample_sizes
                # Series hyperparams
                rec_rate = 0.2
                sd_norm = 0.5
                p_window_size = 30
                linear_basis = (0:0.5:2*pi*cycles)
                n_sample = length(linear_basis)
                series_probs = []
                
                for my_seed in 1:sample_size
                        #Sine wave
                        time_series_sin = map(sin,linear_basis)

                        # White noise  
                        rng = Random.MersenneTwister(my_seed)
                        dist = Distributions.Normal(0,sd_norm)
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
                
                dfs = []
                labels = ["Sine", "WhiteNoise", "SineAndNoise", "AR1low", "AR1hi"]
                
                # "Each plot" level
                for i in 1:5
                        
                        
                        p_vals_all_windows = []
                        avg_vals = []
                        median_vals = []
                        variance_vals = []
                        skewness_vals = []
                        for z in 1:p_window_size
                                p_vals = Float64[]
                                # Get values over samples
                                for j in 1:sample_size
                                        push!(p_vals, series_probs[j][i][z])
                                end
                                
                            # Filter out NaN values
                            p_vals = filter(!isnan, p_vals)
                            #Calculate statistics
                            avg = mean(p_vals)
                            median_val = median(p_vals)
                            variance_val = var(p_vals)
                            skewness_val = skewness(p_vals)
                            # Store values in respective arrays
                            push!(avg_vals, avg)
                            push!(median_vals, median_val)
                            push!(variance_vals, variance_val)
                            push!(skewness_vals, skewness_val)
                            push!(p_vals_all_windows, p_vals)
                        end
                # Create a DataFrame for tabular representation for each i

                df = DataFrame(
                        Label = labels[i],
                        Window = repeat(1:p_window_size, inner=1),
                        Average = avg_vals,
                        Median = median_vals,
                        Variance = variance_vals,
                        Skewness = skewness_vals,
                        )
                push!(dfs,df)
        end
        
        output_df = CSV.write("outputs/Nsamples"*string(sample_size)*"SeriesSize"*string(n_sample)*".csv",vcat(dfs...))
end
end