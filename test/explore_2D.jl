using CSV 
using DataFrames
#using SymbolicInference
#using RecurrenceAnalysis
using Normalization 
using Plots

rec_rate = 0.3;

df = CSV.read("data/sanjuan_data.csv", DataFrame)
ts_raw_p = df[!, "precipitation"]
ts_ma_p = df[!, "mean_precip"]

#  Precipitation

ts_raw_p = convert(Vector{Float64}, ts_raw_p)
ts_ma_p = convert(Vector{Float64}, ts_ma_p[3:end])

norm_fit_p = fit(UnitEnergy, ts_raw_p)
ts_raw_unit_p = normalize(ts_raw_p, norm_fit_p)

norm_max_fit_p = fit(MinMax, ts_raw_p)
ts_raw_max_p = normalize(ts_raw_p, norm_max_fit_p)

norm_fit_ma_p = fit(UnitEnergy, ts_ma_p)
ts_ma_unit_p = normalize(ts_ma_p, norm_fit_ma_p)

norm_max_fit_ma_p = fit(MinMax, ts_ma_p)
ts_ma_max_p = normalize(ts_ma_p, norm_max_fit_ma_p)

# Dengue

ts_raw_d = df[!, "total_cases"]
ts_ma_d = df[!, "mean_cases"]

ts_raw_d = convert(Vector{Float64}, ts_raw_d)
ts_ma_d = convert(Vector{Float64}, ts_ma_d[3:end])

norm_fit_d = fit(UnitEnergy, ts_raw_d)
ts_raw_unit_d = normalize(ts_raw_d, norm_fit_d)

norm_max_fit_d = fit(MinMax, ts_raw_d)
ts_raw_max_d = normalize(ts_raw_d, norm_max_fit_d)

norm_fit_ma_d = fit(UnitEnergy, ts_ma_d)
ts_ma_unit_d = normalize(ts_ma_d, norm_fit_ma_d)

norm_max_fit_ma_d = fit(MinMax, ts_ma_d)
ts_ma_max_d = normalize(ts_ma_d, norm_max_fit_ma_d)


ts_2d = collect(zip(ts_ma_d,ts_ma_p))


Plots.plot(ts_ma_unit_d[1:100],ts_ma_unit_p[1:100],
    arrow=true,marker =:circle, arrowsize=10)

anim = @animate for i in 1:(length(ts_ma_unit_d))
    Plots.plot(ts_ma_unit_d[1:i],ts_ma_unit_p[1:i],
    marker =:circle,xlabel="Dengue cases",
    xlims = (0, 0.4),ylims=(0,0.15),
    ylabel="Precipitation")
end

gif(anim, "anim_fps15.gif", fps = 15)
## Joint Recurrence Matrix
rec_matriz_joint = RecurrenceAnalysis.JointRecurrenceMatrix(
    ts_raw_unit_d, ts_raw_unit_p, rec_rate; fixedrate=true)
rec_matrix_motifs_joint = SymbolicInference.rec_matrix_motifs(
    rec_matriz_joint; max_window=100)

coordinates_joint =  SymbolicInference.extract_recurrences_joint(ts_raw_unit_d, ts_raw_unit_p, rec_matrix_motifs_joint; num_windows=100)

SymbolicInference.plot_motifs_joint(ts_raw_unit_d, ts_raw_unit_p,
    coordinates_joint;n_motifs=3)

#  MA unit
rec_matrix_joint_ma = RecurrenceAnalysis.JointRecurrenceMatrix(ts_ma_unit_d, ts_ma_unit_p, 0.6; fixedrate=true)
rec_matrix_ma_motifs_joint = SymbolicInference.rec_matrix_motifs(rec_matrix_joint_ma; max_window=100)

coordinates_ma_joint =  SymbolicInference.extract_recurrences_joint(ts_ma_unit_d, ts_ma_unit_p, rec_matrix_ma_motifs_joint; num_windows=100)
coordinates_ma_joint
p = SymbolicInference.plot_motifs_joint(ts_ma_unit_d, ts_ma_unit_p,
coordinates_ma_joint[10:end];n_motifs=1)


motif_n = 10
motif_x = vcat(coordinates_ma_joint[motif_n]["x1"],
        coordinates_ma_joint[motif_n]["x2"],
        )

motif_x1 = coordinates_ma_joint[motif_n]["x1"]
motif_x2 = coordinates_ma_joint[motif_n]["x2"]

anim = @animate for i in 1:(length(ts_ma_unit_d))
    Plots.plot(ts_ma_unit_d[1:i],ts_ma_unit_p[1:i],
    marker =:circle,xlabel="Dengue cases",
    xlims = (0, 0.4),ylims=(0,0.15),
    ylabel="Precipitation")

    Plots.plot!(ts_ma_unit_d[filter(x -> x in motif_x1,1:i)],
    ts_ma_unit_p[filter(x -> x in motif_x1,1:i)],
    color="red",linewidth=5)
    Plots.plot!(ts_ma_unit_d[filter(x -> x in motif_x2,1:i)],
    ts_ma_unit_p[filter(x -> x in motif_x2,1:i)],
    color="dark red",linewidth=5)

    Plots.scatter!([ts_ma_unit_d[i]],[ts_ma_unit_p[i]],
    marker =:circle, color="blue")
    if i in motif_x
        Plots.scatter!([ts_ma_unit_d[i]],[ts_ma_unit_p[i]],
        marker =:circle, color="red")
    end

end

gif(anim, "anim3_fps15.gif", fps = 10)
## Cross Recurrence Matrix
rec_matrix_cross = RecurrenceAnalysis.CrossRecurrenceMatrix(ts_raw_max_d, ts_raw_max_p, rec_rate; fixedrate=true)
rec_matrix_motifs_cross = SymbolicInference.rec_matrix_motifs(rec_matrix_cross; max_window=100)

coordinates_cross =  SymbolicInference.extract_recurrences_cross(ts_raw_max_d, ts_raw_max_p, rec_matrix_motifs_cross; num_windows=100)

SymbolicInference.plot_motifs_cross(ts_raw_max_d, ts_raw_max_p,
coordinates_cross;n_motifs=7)


#  MA
rec_matrix_cross_ma = RecurrenceAnalysis.CrossRecurrenceMatrix(ts_ma_max_d, ts_ma_max_p, 0.2; fixedrate=true)
rec_matrix_motifs_cross_ma = SymbolicInference.rec_matrix_motifs(rec_matrix_cross_ma; max_window=100)

coordinates_cross_ma =  SymbolicInference.extract_recurrences_cross(ts_ma_max_d, ts_ma_max_p, 
rec_matrix_motifs_cross_ma; num_windows=100)

SymbolicInference.plot_motifs_cross(ts_ma_max_d, ts_ma_max_p,
coordinates_cross_ma;n_motifs=2)