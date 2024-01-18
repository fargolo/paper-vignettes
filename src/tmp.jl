using TaylorSeries

function p_val_weighted(p,q,l,n)

    if (p+q - 1 > 0.01) || (l > n)
        println("p + q must be equal to 1 and l <= n")
            return(NaN)
    end 

    weighted_coeffs = map(x->weighted_bin_runs_coeff(p,q,x,n),(l-1):1:n)
    probs = diff(weighted_coeffs)
    Float64(sum(probs))

end

function weighted_bin_runs_coeff(p,q,l,n)

    if (p+q - 1 > 0.01) || (l > n)
        println("p + q must be equal to 1 and l <= n")
            return(NaN)
    end 

    wei_mgf(p,q,l,z) = (1 - p^l * z^l )/(1 - z + q*(p^l)*(z^(l+1)))
    tay_exp = TaylorSeries.taylor_expand(z -> wei_mgf(p,q,l,z),order=n+1)
    TaylorSeries.getcoeff(tay_exp,n)

end
