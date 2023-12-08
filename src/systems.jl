
"""

Random Gaussian (N(0,1)) walk; c scales increments.
x_0 = x ~ N(0,1)
x_t = x_{t-1} + c*e , e ~ N(0,1) 
"""
function rand_walk_gauss(ts_size;dist = Normal(0,1),c=1)
        ts = []
        push!(ts,rand(dist,1)[1])
        for i in 2:ts_size
                cur_val = ts[i-1] + c*rand(dist,1)[1]
                push!(ts,cur_val)
        end
        Float64.(ts)
end

