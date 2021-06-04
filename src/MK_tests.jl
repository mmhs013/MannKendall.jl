function mk_original_test(x, α = 0.05)
    n = length(x)
    s = mk_score(x)
    var_s = variance_s(x,length(x))
    τ = s/(.5*n*(n-1))
    z = z_score(s, var_s)
    
    p, h, trend = p_value(z, α)
    
    return (trend=trend,  h=h, p=p, z=z, τ=τ, s=s, var_s=var_s)
end