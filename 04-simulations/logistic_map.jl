using SymbolicInference
using Random, Distributions
using RecurrenceAnalysis
using Statistics , StatsBase
using DataFrames , CSV
using DynamicalSystems
using PredefinedDynamicalSystems


# 4.669201609 = (Rn-Rn-1)/(Rn+1-Rn)

logistic_ds2 = PredefinedDynamicalSystems.logistic(0.1; r = 3.2)
logistic_ds4 = PredefinedDynamicalSystems.logistic(0.1; r = 3.5)
logistic_ds8 = PredefinedDynamicalSystems.logistic(0.1; r = 3.56)
logistic_ds_chaos = PredefinedDynamicalSystems.logistic(0.1; r = 3.5699456)

total_time = 60
X2, t2 = trajectory(logistic_ds2, total_time)
X4, t4 = trajectory(logistic_ds4, total_time)

Plots.plot([Matrix(X2)...])
Plots.plot([Matrix(X4)...])


SymbolicInference.persistence_barcode([Matrix(X2)...]; 
        range = collect(0.3:0.3:0.9),
        n_windows=1,alpha_thresh=0.01)


rm_r32_thr005 = RecurrenceAnalysis.RecurrenceMatrix(X2,0.05;
                fixedrate=true)
rm_r35_thr005 = RecurrenceAnalysis.RecurrenceMatrix(X4,0.05;
                fixedrate=true)

rm_r32_thr005_mot = rec_matrix_motifs(rm_r32_thr005;
        window_range=[2,4,6,8,10,12,14,16],n_motifs=1)
rm_r35_thr005_mot = rec_matrix_motifs(rm_r35_thr005;
        window_range=[4,8,12,16,20],n_motifs=1)


coords_rm_r32_thr005 = SymbolicInference.extract_recurrences([Matrix(X2)...], 
                rm_r32_thr005_mot; num_windows=8)

coords_rm_r35_thr005 = SymbolicInference.extract_recurrences([Matrix(X4)...], 
                rm_r35_thr005_mot; num_windows=5)

#using GLMakie
#GLMakie.activate!()
## 7,8 ; 19,20 ; 21,22 / 150,151 / 219 valley
p2 = SymbolicInference.plot_motifs([Matrix(X2)...], coords_rm_r32_thr005; n_motifs=4)
p4 = SymbolicInference.plot_motifs([Matrix(X4)...], coords_rm_r35_thr005[4:5]; n_motifs=1)

r = 3.57 # Transition to chaos

include("systems.jl")

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

                       all_probs = map(x-> double_inference_weighted(x;
                                        max_window=p_window_size,seqs="recurrences"),res_recs)
                    
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

# using JLD2
# jldsave("outputs/full_sim.jld2"; series_probs)