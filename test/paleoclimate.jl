using Pkg
Pkg.add(path="/home/boitata/repos/SymbolicInference.jl")
using SymbolicInference
using Plots

include("/home/boitata/repos/paper-vignettes/src/VignetteDependencies.jl")

rec_rate=0.2

folder_path = "data/paleoclimate/short_series"
outputs_path = "outputs/paleoclimate/short_series"
file_list = filter(x -> endswith(x, ".txt"), readdir(folder_path))


# Short

series_probs = []
for file_name in file_list
    file_path = joinpath(folder_path, file_name)
    file_output_path = joinpath(outputs_path,file_name*"double.png") 
    
    df = CSV.read(file_path,DataFrame;
                delim=" ", header=["Index","Time","Data"], skipto=2)
    
    data_array = df[:, end]
    
    println("Data from $file_name:")
    println(data_array)


    res_rec = RecurrenceAnalysis.RecurrenceMatrix(data_array,rec_rate;fixedrate=true)
    
    p_window_size = length(data_array) - 3

    all_probs = double_inference_weighted(res_rec;max_window=p_window_size,seqs="double")
    
    p_plt = scatter(1:length(all_probs),all_probs,ylims=(0,1))

    savefig(p_plt,file_output_path)

    push!(series_probs,all_probs)

end


# Medium

rec_rate=0.2

folder_path = "data/paleoclimate/medium_series"
outputs_path = "outputs/paleoclimate/medium_series"
file_list = filter(x -> endswith(x, ".txt"), readdir(folder_path))


series_probs = []
for file_name in file_list
    file_path = joinpath(folder_path, file_name)
    file_output_path = joinpath(outputs_path,file_name*".png") 
    
    df = CSV.read(file_path,DataFrame;
                delim=" ", header=["Index","Time","Data"], skipto=2)
    
    data_array = df[:, end]
    
    println("Data from $file_name:")
    println(data_array)


    res_rec = RecurrenceAnalysis.RecurrenceMatrix(data_array,rec_rate;fixedrate=true)
    
    p_window_size = length(data_array) - 3

    all_probs = double_inference_weighted(res_rec;max_window=p_window_size,seqs="recurrences")
    
    p_plt = scatter(1:length(all_probs),all_probs,ylims=(0,1))

    savefig(p_plt,file_output_path)

    push!(series_probs,all_probs)

end


# Long

rec_rate=0.2

folder_path = "data/paleoclimate/long_series"
outputs_path = "outputs/paleoclimate/long_series"
file_list = filter(x -> endswith(x, ".txt"), readdir(folder_path))


series_probs = []
for file_name in file_list
    file_path = joinpath(folder_path, file_name)
    file_output_path = joinpath(outputs_path,file_name*".png") 
    
    df = CSV.read(file_path,DataFrame;
                delim=" ", header=["Index","Time","Data"], skipto=2)
    
    data_array = df[:, end]
    
    println("Data from $file_name:")
    println(data_array)


    res_rec = RecurrenceAnalysis.RecurrenceMatrix(data_array,rec_rate;fixedrate=true)
    
    p_window_size = length(data_array) - 3

    all_probs = double_inference_weighted(res_rec;max_window=p_window_size,seqs="recurrences")
    
    p_plt = scatter(1:length(all_probs),all_probs,ylims=(0,1))

    savefig(p_plt,file_output_path)

    push!(series_probs,all_probs)

end
