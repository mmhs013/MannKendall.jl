function mk_score(x)
    s = 0
    n = length(x)
    
    for i = 1:n-1
        for j = i : n
            s = s + sign(x[j] - x[i])
        end
    end
    
    return s
end

function variance_s(x, n)
    # calculate the unique data
    unique_x = unique(x)
    g = length(unique_x)

    # calculate the var(s)
    if n == g            # there is no tie
        var_s = (n*(n-1)*(2*n+5))/18
        
    else                 # there are some ties in data
        tp = zeros(length(unique_x))
        demo = ones(n)
        
        for i in 1:g
            tp[i] = sum(demo[x .== unique_x[i]])
        end
        var_s = (n*(n-1)*(2*n+5) - sum(tp*(tp-1)*(2*tp+5)))/18
    end
        
    return var_s
end

function z_score(s, var_s)
    if s > 0
        z = (s - 1)/sqrt(var_s)
    elseif s == 0
        z = 0
    elseif s < 0
        z = (s + 1)/sqrt(var_s)
    end
    
    return z
end

function p_value(z, α)
    # two tail test
    p = 2*(1-normcdf(abs(z)))  
    h = abs(z) > norminvcdf(1-α/2)

    if (z < 0) && h
        trend = "decreasing"
    elseif (z > 0) && h
        trend = "increasing"
    else
        trend = "no trend"
    end
    
    return p, h, trend
end
