using Plots
using StatsPlots

# Plot recurrence matrices

# Series

## Hyperparams
cycles = 8 # cycles on sin wave
sd_norm = 0.5
p_window_size = 50

## Sine wave 
time_series_sin = map(sin,0:0.5:2*pi*cycles)

## White noise  
rng = Random.MersenneTwister(my_seed)
n_sample=length(time_series_sin)
dist = Distributions.Normal(0,sd_norm)
white_n_time_series = rand(dist,n_sample)

## Noisy sine wave
time_series_sin_nois = white_n_time_series .+ time_series_sin

## Gaussian random walk (AR1)
## https://discourse.julialang.org/t/simple-random-walk-performance-tips/61553/2
rand_walk_gauss_ts50 = rand_walk_gauss(n_sample;a=0.5,c=0.5)
rand_walk_gauss_ts90 = rand_walk_gauss(n_sample;a=0.9,c=0.5)

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


# Recurrence matrices 

res_recs = map(x -> RecurrenceAnalysis.RecurrenceMatrix(x,0.2;fixedrate=true),
[time_series_sin,white_n_time_series,time_series_sin_nois,rand_walk_gauss_ts50,rand_walk_gauss_ts90])
res_recs2 = map(x -> RecurrenceAnalysis.RecurrenceMatrix(x,0.5;fixedrate=true),
        [time_series_sin,white_n_time_series,time_series_sin_nois,rand_walk_gauss_ts50,rand_walk_gauss_ts90])
res_recs3 = map(x -> RecurrenceAnalysis.RecurrenceMatrix(x,0.8;fixedrate=true),
        [time_series_sin,white_n_time_series,time_series_sin_nois,rand_walk_gauss_ts50,rand_walk_gauss_ts90])


myscheme = ColorScheme([Colors.RGB(0.94118, 0.60784, 0.58823), 
                        Colors.RGB(0.67058, 0.3843137, 0.3607843)],
               "custom", "twotone, red and green")



loadcolorscheme(:myscheme,ColorScheme([get(myscheme, i) for i in 0.0:0.01:1.0]))
ps = map(x -> heatmap(x,legend = false,c= :myscheme),res_recs)
p_recs = plot(ps...;size= (2560,512),layout=(1,5))
ps2 = map(x -> heatmap(x,legend = false,c= :starrynight),res_recs2)
p_recs2 = plot(ps2...;size= (2560,512),layout=(1,5))
ps3 = map(x -> heatmap(x,legend = false,c= :starrynight),res_recs3)
p_recs3 = plot(ps3...;size= (2560,512),layout=(1,5))


p_series_recs = plot(p_series,p_recs,p_recs2,p_recs3,layout=(4,1),size= (2560,2048))


# Boxplot 

## Get probabilities
# Use code in "generate_probs.jl" to generate a series_probs object

p_plots = []
for i in 1:5
    p = plot()
    # "Each window" level
    p_vals_all_windows = []
    for z in 1:30
        p_vals = Float64[]
        # Get values over samples
        for j in 1:1000
            push!(p_vals,series_probs[j][i][z])
        end
        println("p_vals is: ",p_vals)
        push!(p_vals_all_windows,filter(!isnan,(p_vals)))
    end
    try
        boxplot!(p_vals_all_windows,size=(768,256),fa=0.5;outliers=false)
        #dotplot!(p_vals_all_windows,size=(768,256))
        #violin!(p_vals_all_windows,size=(768,256))
   
    catch e
        println("Nan")
    end
    push!(p_plots,p)
end

p_boxs = plot(p_plots...,layout=(1,5),size=(2048,256),legend=false)

l = @layout [
     a 
     b{0.2h}]
             
plot(p_series_recs,p_boxs,layout=l)