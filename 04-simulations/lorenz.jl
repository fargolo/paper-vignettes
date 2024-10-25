using SymbolicInference
using Random, Distributions
using RecurrenceAnalysis
using Statistics , StatsBase
using DataFrames , CSV
using DynamicalSystems
using PredefinedDynamicalSystems


# 4.669201609 = (Rn-Rn-1)/(Rn+1-Rn)

lorenz_ds = PredefinedDynamicalSystems.lorenz([0.0, 10.0, 0.0]; σ = 10.0, ρ = 28.0, β = 8/3)
total_time = 300
X, t = trajectory(lorenz_ds, total_time)

plot(Matrix(X)[:,1], 
        Matrix(X)[:,2],
        Matrix(X)[:,3])


#SymbolicInference.persistence_barcode([Matrix(X2)...]; 
#        range = collect(0.3:0.3:0.9),
#        n_windows=1,alpha_thresh=0.01)


lorenz_rm = RecurrenceAnalysis.RecurrenceMatrix(Matrix(X),
                        0.05;fixedrate=true)



rm_r32_thr005_mot = rec_matrix_motifs(rm_r32_thr005;
        window_range=collect(1:40),n_motifs=1)
rm_r35_thr005_mot = rec_matrix_motifs(rm_r35_thr005;
        window_range=collect(1:40),n_motifs=1)
rm_r356_thr005_mot = rec_matrix_motifs(rm_r356_thr0056;
        window_range=collect(1:40),n_motifs=1)

bin_miss(x) = ifelse(ismissing(x),true,false)

a = plot(bin_miss.(rm_r32_thr005_mot["Motifs starts and duration"]))
b = plot(bin_miss.(rm_r35_thr005_mot["Motifs starts and duration"]))
c = plot(bin_miss.(rm_r356_thr005_mot["Motifs starts and duration"]))

logis_plot = plot(a,b,c,layout = (3, 1),legend=false,
        title=["r=3.2" "r=3.5" "r = 3.5699456"])

#savefig(logis_plot,"logistic_plot.png")

coords_rm_r32_thr005 = SymbolicInference.extract_recurrences([Matrix(X2)...], 
                rm_r32_thr005_mot; num_windows=8)

coords_rm_r35_thr005 = SymbolicInference.extract_recurrences([Matrix(X4)...], 
                rm_r35_thr005_mot; num_windows=5)

#using GLMakie
#GLMakie.activate!()
## 7,8 ; 19,20 ; 21,22 / 150,151 / 219 valley
p2 = SymbolicInference.plot_motifs([Matrix(X2)...], coords_rm_r32_thr005; n_motifs=4)
p4 = SymbolicInference.plot_motifs([Matrix(X4)...], coords_rm_r35_thr005[4:5]; n_motifs=1)

