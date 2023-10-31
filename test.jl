using Plots

include("src/SymbolicInference.jl")

# Hyperparams
my_seed = 1234
cycles = 8
sd_norm = 1
rec_rate = 0.5

rng = MersenneTwister(my_seed)

time_series_sin = map(sin,0:0.3:2*pi*cycles)
n_sample=length(time_series_sin)
dist = Normal(0,sd_norm)
time_series = rand(dist,n_sample)
time_series_sin_nois = time_series .+ time_series_sin
p_tsn = plot(1:length(time_series_sin_nois),time_series_sin_nois)
p_sin = plot(1:length(time_series_sin_nois),time_series_sin)
p_r = plot(1:length(time_series_sin_nois),time_series)
plot(p_tsn,p_sin,p_r,layout=(3,1))
res_recs = map(x -> RecurrenceAnalysis.RecurrenceMatrix(x,rec_rate;fixedrate=true),
        [time_series_sin_nois,time_series_sin,time_series])
        res_infs = map(x -> double_inference(x,rec_rate;fixedrate=true),res_recs)


time_series_sin_nois = RecurrenceAnalysis.RecurrenceMatrix(time_series_sin_nois,rec_rate;fixedrate=true)
time_series_sin = RecurrenceAnalysis.RecurrenceMatrix(time_series_sin,rec_rate;fixedrate=true)
time_series = RecurrenceAnalysis.RecurrenceMatrix(time_series,rec_rate;fixedrate=true)


results_sin_r = map(x -> double_inference(time_series_sin_nois;seqs=x),["double","poincare","recurrences"])
results_sin = map(x -> double_inference(time_series_sin;seqs=x),["double","poincare","recurrences"])
results_r = map(x -> double_inference(time_series;seqs=x),["double","poincare","recurrences"])


ps_sin_r = map(plot,results_sin_r)
ps_sin = map(plot,results_sin)
ps_r = map(plot,results_r)

p_sin_r = plot(ps_sin_r[1],ps_sin_r[2],ps_sin_r[3],layout=(3,1);ylims = (0, 0.3),title="sin(x) + N(0,3)")
p_sin = plot(ps_sin[1],ps_sin[2],ps_sin[3],layout=(3,1);ylims = (0, 0.3),title="sin(x)")
p_r = plot(ps_r[1],ps_r[2],ps_r[3],layout=(3,1);ylims = (0, 0.3),title="N(0,3)")

plot(p_sin_r,p_sin,p_r,layout=(1,3))
plot!(size=(1600,1600))
plot(p_sin_r[3],p_sin[3],p_r[3]; title = "Sin(x) vs. N(0,1)")

