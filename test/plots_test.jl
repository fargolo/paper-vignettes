using Plots

include("test_v1.jl")
# Assuming series_probs is your vector of vectors of vectors
# Replace this with your actual data
#series_probs = [randn(50) for _ in 1:5, _ in 1:100]

# Create a plot for each element in the second dimension

# Plot level (e.g. sine, noise, sine noise...)
p_plots = []
for i in 1:5
    p = plot()
    # "Each window" level
    p_vals_all_windows = []
    for z in 1:50
        p_vals = Float64[]
        # Get values over samples
        for j in 1:300
            push!(p_vals,series_probs[j][i][z])
        end
        push!(p_vals_all_windows,filter(!isnan,(p_vals)))
    end
    try
        boxplot!(p_vals_all_windows,fa=0.5)
        
    catch e
        println("Nan")
    end
    push!(p_plots,p)
end

plot(p_plots...,layout=(5,1),size=(2560,512),legend=false)

plot(p_plots[3],legend=false)

skipmissing(series_probs[1][1])

for i in 1:length(series_probs[1][1])
    plot()

    # Iterate over the first dimension (size 100)
    for j in 1:length(series_probs)
        # Extract data for the current boxplot
        data = hcat([series_probs[j][k][i] for k in 1:length(series_probs[j]) if !isnan(series_probs[j][k][i])]...)

        # Add a boxplot to the current plot
        boxplot!(data, label="Boxplot $j")
    end

    # Customize plot labels and title
    xlabel!("Boxplots")
    ylabel!("Values")
    title!("Boxplots for Series $i")
end

# Display the plots
display(Plots.plot!())