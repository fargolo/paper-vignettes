using Random
using AnalyticComb
using StatsBase
using Plots

n_series , n_sim = 200 , 1000 
probs = []
for i in 1:n_sim
    cur_array = Random.bitrand(n_series)
    arr_enc = StatsBase.rle(cur_array)
    max_val = maximum(arr_enc[2])
    cur_prob = AnalyticComb.weighted_bin_runs_pval(0.5,0.5,max_val,n_series)
    push!(probs,cur_prob)
end

alpha_thresh = collect(0:0.01:1)
props_p = []
for cur_thresh in alpha_thresh
    pos_vals = map(x -> ifelse(x<cur_thresh,1,0),probs)
    alpha_prop_sim = sum(pos_vals)/length(pos_vals)
    push!(props_p,alpha_prop_sim)
end

scatter(1:length(props_p),props_p)
scat_pl = scatter(1:length(props_p),props_p,
    xticks=(collect(1:2:length(alpha_thresh)),
    string.(alpha_thresh)[1:2:length(alpha_thresh)]),
    xrotation=90)
plot!(collect(1:length(alpha_thresh)),
    collect(0:0.01:1))
ecdf_func = StatsBase.ecdf(convert.(Float64,props_p))
plot(ecdf_func(0:0.025:1))

edcf = StatsBase.ecdf(convert.(Float64,probs))
edcf

lin_pl = Plots.plot(sort(probs))
lin_pl2 = Plots.plot(sort(probs), (1:1000)./1000)
ecdf_func = StatsBase.ecdf(convert.(Float64,probs))
plot(ecdf_func(0:0.02:1))

pos_vals = map(x -> ifelse(x<0.05,1,0),probs)
alpha_prop_sim = sum(pos_vals)/length(pos_vals)

Plots.savefig(scat_pl,"outputs/line_plot_p-values4.png")