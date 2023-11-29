using Distributed
using SharedArrays
using Plots
include("src/SymbolicInference.jl") 
#addprocs(2)

my_seed = 1234
cycles = 6
sd_norm = 0.5
rec_rate = 0.2

rng = MersenneTwister(my_seed)

time_series_sin = map(sin,0:0.3:2*pi*cycles)

n_sample=length(time_series_sin)
dist = Normal(0,sd_norm)
time_series = rand(dist,n_sample)

time_series_sin_nois = time_series .+ time_series_sin


p_tsn = plot(1:length(time_series_sin_nois),time_series_sin_nois)
p_sin = plot(1:length(time_series_sin_nois),time_series_sin)
p_r = plot(1:length(time_series_sin_nois),time_series)
p_3_series = plot(p_tsn,p_sin,p_r,layout=(3,1))


@everywhere begin
        using Pkg;
        include("src/SymbolicInference.jl") 
        Pkg.precompile()
end
    
@everywhere begin 
        using Plots

        my_seed = 1234
        cycles = 5
        sd_norm = 0.5
        rec_rate = 0.2

        rng = MersenneTwister(my_seed)

        time_series_sin = map(sin,0:0.3:2*pi*cycles)
        n_sample=length(time_series_sin)
        dist = Normal(0,sd_norm)
        time_series = rand(dist,n_sample)
        time_series_sin_nois = time_series .+ time_series_sin
end

@everywhere res_recs = pmap(x -> RecurrenceAnalysis.RecurrenceMatrix(x,rec_rate;fixedrate=true),
                [time_series_sin_nois,time_series_sin,time_series])


#res_infs = SharedArray[]
#@everywhere for i in ["double","poincare","recurrence"]
#        cur_res = pmap(x -> double_inference(x;seqs = i),res_recs)
#        push!(res_infs,cur_res)
#end

res_infs_d = pmap(x -> double_inference_weighted(x;seqs = "double"),res_recs)
res_infs_p = pmap(x -> double_inference_weighted(x;seqs = "poincare"),res_recs)
res_infs_r = pmap(x -> double_inference_weighted(x;seqs = "recurrences"),res_recs)

ps_sin_r = map(plot,res_infs_d)
ps_sin = map(plot,res_infs_p)
ps_r = map(plot,res_infs_r)

p_sin_r = plot(ps_sin_r[1],ps_sin_r[2],ps_sin_r[3],layout=(3,1);ylims = (0, 0.3),title="Double")
p_sin = plot(ps_sin[1],ps_sin[2],ps_sin[3],layout=(3,1);ylims = (0, 0.3),title="Poincare")
p_r = plot(ps_r[1],ps_r[2],ps_r[3],layout=(3,1);ylims = (0, 0.3),title="Recurrence")

p_3_hypothesis = plot(p_sin_r,p_sin,p_r,layout=(1,3),size=(1600,1600))

l = @layout [
    a{0.3w} b{0.7w}
]

plot(p_3_series,p_3_hypothesis,layout=l)
