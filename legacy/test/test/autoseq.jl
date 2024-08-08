
rec_rate = 0.2
p_window_size = 38


function parse_numbers(s)
    pieces = split(s, ' ', keepempty=false)
    map(pieces) do piece
        parse(Float64, piece)
    end
end

# Automatic sequence

folder_path = "data/autoseq"
file_name = "seq5.txt"
file_path = joinpath(folder_path, file_name)

str_sq = read(file_path,String)


seq_series = parse_numbers(str_sq)

res_rec = RecurrenceAnalysis.RecurrenceMatrix(seq_series,rec_rate;fixedrate=true)
all_probs = double_inference_weighted(res_rec;max_window=p_window_size,seqs="recurrences")

heatmap(res_rec)
heatmap(res_rec[709:729,1:20])
heatmap(res_rec[1:20,1:20])

p_plt = scatter(1:length(all_probs),all_probs,ylims=(0,1))

# Fib

folder_path = "data/autoseq"
file_name = "fib.txt"
file_path = joinpath(folder_path, file_name)

str_sq = readline(file_path)

seq_series = split(str_sq,"")
seq_series = parse.(Int64,seq_series)

res_rec = RecurrenceAnalysis.RecurrenceMatrix(seq_series,rec_rate;fixedrate=true)
all_probs = double_inference_weighted(res_rec;max_window=p_window_size,seqs="recurrences")

heatmap(res_rec)
heatmap(res_rec[709:729,1:20])
heatmap(res_rec[1:20,1:20])

p_plt = scatter(1:length(all_probs),all_probs,ylims=(0,1))


# Period

folder_path = "data/autoseq"
file_name = "period.txt"
file_path = joinpath(folder_path, file_name)

str_sq = readline(file_path)

seq_series = split(str_sq,"")
seq_series = parse.(Int64,seq_series)

res_rec = RecurrenceAnalysis.RecurrenceMatrix(seq_series,rec_rate;fixedrate=true)
all_probs = double_inference_weighted(res_rec;max_window=p_window_size,seqs="recurrences")

heatmap(res_rec)
heatmap(res_rec[709:729,1:20])
heatmap(res_rec[1:20,1:20])

p_plt = scatter(1:length(all_probs),all_probs,ylims=(0,1))
