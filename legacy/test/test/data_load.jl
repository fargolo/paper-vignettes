using CSV 
using DataFrames
using Normalization 

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
