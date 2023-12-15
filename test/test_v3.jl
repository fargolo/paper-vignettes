using Plots

include("../src/SymbolicInference.jl")

# Time-series types:  
## Sine wave
## White noise
## Noisy sine wave
## Gaussian random walk (AR1)

# Hyperparams
my_seed = 1234
cycles = 8 # cycles on sin wave
rec_rate = 0.5
p_window_size = 10

# Sine wave 

time_series_sin = map(sin,0:0.7:2*pi*cycles)

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

## Sine wave
p_sin = plot(1:length(time_series_sin_nois),time_series_sin,title="Sine wave")
## White noise
p_wn = plot(1:length(time_series_sin_nois),white_n_time_series,title="White noise")
## Noisy sine wave
p_tsn = plot(1:length(time_series_sin_nois),time_series_sin_nois,title="Sine wave + white noise")
## Gaussian random walk (AR1)
p_rwg50 = plot(1:length(time_series_sin_nois),rand_walk_gauss_ts50,title="Gauss Rand walk - N(0,1) ; a=0.5,c=0.5")
p_rwg90 = plot(1:length(time_series_sin_nois),rand_walk_gauss_ts90,title="Gauss Rand walk - N(0,1) ; a=0.9,c=0.5")

p_series = plot(p_sin,p_wn,p_tsn,p_rwg50,p_rwg90,layout=(1,5),
        size= (2560,512))

res_recs = map(x -> RecurrenceAnalysis.RecurrenceMatrix(x,rec_rate;fixedrate=true),
        [time_series_sin,white_n_time_series,time_series_sin_nois,rand_walk_gauss_ts50,rand_walk_gauss_ts90])

ps = map(x -> heatmap(x,legend = false,c= :starrynight),res_recs)
p_recs = plot(ps...;size= (2560,512),layout=(1,5))

plot(p_series,p_recs,layout=(2,1))

res_recs = map(x -> RecurrenceAnalysis.RecurrenceMatrix(x,0.2;fixedrate=true),
        [time_series_sin,white_n_time_series,time_series_sin_nois,rand_walk_gauss_ts50,rand_walk_gauss_ts90])

ps = map(x -> heatmap(x,legend = false,c= :starrynight),res_recs)
p_recs2 = plot(ps...;size= (2560,512),layout=(1,5))

res_recs = map(x -> RecurrenceAnalysis.RecurrenceMatrix(x,0.8;fixedrate=true),
        [time_series_sin,white_n_time_series,time_series_sin_nois,rand_walk_gauss_ts50,rand_walk_gauss_ts90])

ps = map(x -> heatmap(x,legend = false,c= :starrynight),res_recs)
p_recs3 = plot(ps...;size= (2560,512),layout=(1,5))


plot(p_series,p_recs2,p_recs,p_recs3,layout=(4,1),
size= (2560,2048))


rqa(res_recs[1])

all_probs = map(x-> double_inference_weighted(x;max_window=p_window_size,seqs="recurrences"),res_recs)

p_vals_plot = map(x -> plot(collect(1:p_window_size),x,seriestype=:scatter),all_probs)
p_vals_plot = plot(p_vals_plot...;size= (2560,512),layout=(1,5))
ylims!(0,1)

plot(p_series,p_recs2,p_recs,p_recs3,p_vals_plot,layout=(5,1),
size= (2560,2048))


#n = size(res_recs[1])[1]
#[ res_recs[1][i,i] for i=1:n ]
#[ res_recs[1][n-i+1,i] for i=1:n ]
#res_infs = map(x -> double_inference_weighted(x),res_recs)


#time_series_sin_nois = RecurrenceAnalysis.RecurrenceMatrix(time_series_sin_nois,rec_rate;fixedrate=true)
#time_series_sin = RecurrenceAnalysis.RecurrenceMatrix(time_series_sin,rec_rate;fixedrate=true)
#time_series = RecurrenceAnalysis.RecurrenceMatrix(time_series,rec_rate;fixedrate=true)


#results_sin_r = map(x -> double_inference(time_series_sin_nois;seqs=x),["double","poincare","recurrences"])
#results_sin = map(x -> double_inference(time_series_sin;seqs=x),["double","poincare","recurrences"])
#results_r = map(x -> double_inference(time_series;seqs=x),["double","poincare","recurrences"])


#ps_sin_r = map(plot,results_sin_r)
#ps_sin = map(plot,results_sin)
#ps_r = map(plot,results_r)

#p_sin_r = plot(ps_sin_r[1],ps_sin_r[2],ps_sin_r[3],layout=(3,1);ylims = (0, 0.3),title="sin(x) + N(0,3)")
#p_sin = plot(ps_sin[1],ps_sin[2],ps_sin[3],layout=(3,1);ylims = (0, 0.3),title="sin(x)")
#p_r = plot(ps_r[1],ps_r[2],ps_r[3],layout=(3,1);ylims = (0, 0.3),title="N(0,3)")

#plot(p_sin_r,p_sin,p_r,layout=(1,3))
#plot!(size=(1600,1600))
#plot(p_sin_r[3],p_sin[3],p_r[3]; title = "Sin(x) vs. N(0,1)")

