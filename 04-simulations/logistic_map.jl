using SymbolicInference
using Plots
using Random, Distributions
using RecurrenceAnalysis
using Statistics , StatsBase
using DataFrames , CSV
using DynamicalSystems
using PredefinedDynamicalSystems


# 4.669201609 = (Rn-Rn-1)/(Rn+1-Rn)

logistic_ds2 = PredefinedDynamicalSystems.logistic(0.1; r = 3.2)
logistic_ds4 = PredefinedDynamicalSystems.logistic(0.1; r = 3.5)
logistic_ds_trans = PredefinedDynamicalSystems.logistic(0.1; r = 3.56)
logistic_ds_chaos = PredefinedDynamicalSystems.logistic(0.1; r = 4)

total_time = 300
X2, t2 = trajectory(logistic_ds2, total_time)
X4, t4 = trajectory(logistic_ds4, total_time)
Xt, tt = trajectory(logistic_ds_trans, total_time)
Xc, tc = trajectory(logistic_ds_chaos, total_time)


a1 = Plots.plot([Matrix(X2)...][1:200])
b1 = Plots.plot([Matrix(X4)...][1:200])
c1 = Plots.plot([Matrix(Xt)...][1:200])
d1 = Plots.plot([Matrix(Xc)...][1:200])
traj_plot = plot(a1,b1,c1,d1, layout = (4, 1),legend=false,
        title=["r = 3.20" "r = 3.50" "r = 3.56" "r = 4.00"])



SymbolicInference.persistence_barcode([Matrix(Xc)...]; 
        range = collect(0.1:0.1:0.9), n_windows=1, alpha_thresh=0.05)

rm_r32_thr005 = RecurrenceAnalysis.RecurrenceMatrix(X2,0.05;fixedrate=true)

rm_r35_thr005 = RecurrenceAnalysis.RecurrenceMatrix(X4,0.05;fixedrate=true)

rm_r357_thr005 = RecurrenceAnalysis.RecurrenceMatrix(Xt,0.05;fixedrate=true)

rm_r4_thr0056 = RecurrenceAnalysis.RecurrenceMatrix(Xc,0.05; fixedrate=true)



rm_r32_thr005_mot = rec_matrix_motifs(rm_r32_thr005; window_range=collect(1:50),n_motifs=1)

rm_r35_thr005_mot = rec_matrix_motifs(rm_r35_thr005; window_range=collect(1:50),n_motifs=1)

rm_r357_thr005_mot = rec_matrix_motifs(rm_r357_thr005; window_range=collect(1:50),n_motifs=1)

rm_r4_thr0056_mot = rec_matrix_motifs(rm_r4_thr0056; window_range=collect(1:50),n_motifs=1)

bin_miss(x) = ifelse(ismissing(x),false,true)

a = plot(bin_miss.(rm_r32_thr005_mot["Motifs starts and duration"]))
b = plot(bin_miss.(rm_r35_thr005_mot["Motifs starts and duration"]))
c = plot(bin_miss.(rm_r357_thr005_mot["Motifs starts and duration"]))
d = plot(bin_miss.(rm_r4_thr0056_mot["Motifs starts and duration"]),
                ylims = (0,1))


logis_plot = plot(a,b,c,d, layout = (4, 1),legend=false,
        title=["r = 3.20" "r = 3.50" "r = 3.56" "r = 4.00"])

multiplot = plot(traj_plot,logis_plot,layout = (1, 2))

#savefig(multiplot,"logistic_plot.png")

coords_rm_r32_thr005 = SymbolicInference.extract_recurrences([Matrix(X2)...], 
                rm_r32_thr005_mot; num_windows=8)

coords_rm_r35_thr005 = SymbolicInference.extract_recurrences([Matrix(X4)...], 
                rm_r35_thr005_mot; num_windows=5)

#using GLMakie
#GLMakie.activate!()
## 7,8 ; 19,20 ; 21,22 / 150,151 / 219 valley
p2 = SymbolicInference.plot_motifs([Matrix(X2)...], coords_rm_r32_thr005; n_motifs=4)
p4 = SymbolicInference.plot_motifs([Matrix(X4)...], coords_rm_r35_thr005[4:5]; n_motifs=1)

