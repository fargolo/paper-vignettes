using Plots

include("src/SymbolicInference.jl")

# Hyperparams
my_seed = 1234
cycles = 8
sd_norm = 3
rec_rate = 0.5

rng = MersenneTwister(my_seed)

time_series_sin = map(sin,0:0.3:2*pi*cycles)
n_sample=length(time_series_sin)
dist = Normal(0,sd_norm)
time_series = rand(dist,n_sample)

time_series = RecurrenceAnalysis.RecurrenceMatrix(time_series,rec_rate;fixedrate=true)
time_series_sin = RecurrenceAnalysis.RecurrenceMatrix(time_series_sin,rec_rate;fixedrate=true)

results = map(x -> double_inference(time_series_sin;seqs=x),["double","poincare","recurrences"])
results_r = map(x -> double_inference(time_series;seqs=x),["double","poincare","recurrences"])

ps = map(plot,results)
ps_r = map(plot,results_r)

p_sin = plot(ps[1],ps[2],ps[3],layout=(3,1);ylims = (0, 0.3),title="sin(x)")
p_r = plot(ps_r[1],ps_r[2],ps_r[3],layout=(3,1);ylims = (0, 0.3),title="N(0,1)")

plot(p_sin,p_r,layout=(1,2))
plot(ps[3],ps_r[3]; title = "Sin(x) vs. N(0,1)")

