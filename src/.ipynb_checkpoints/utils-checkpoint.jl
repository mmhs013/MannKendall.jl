"""
    __preprocessing(x_old)

Internal preprocessing routine shared by the Mann-Kendall tests and Sen's slope
estimators. It mirrors ``__preprocessing`` from pyMannKendall.

The input is converted to a ``Float64`` container.

Returns a tuple ``(x, c)`` where:

- ``x`` is the data as a ``Vector{Float64}`` for univariate input, or as a
  ``Matrix{Float64}`` with ``c`` columns for multivariate input (when the matrix
  has a single column it is flattened to a vector, exactly like pyMannKendall).
- ``c`` is the number of columns / series.

Only explicit ``for`` loops are used here (no vectorised broadcasting) because
plain loops are fast and idiomatic in Julia.
"""
function __preprocessing(x_old)
    if x_old isa AbstractMatrix
        x = Float64.(x_old)
        (n, c) = size(x)
        if c == 1
            x = vec(x)
        end
        return (x, c)
    else
        x = Float64.(collect(x_old))
        return (x, 1)
    end
end

"""
    __missing_values_analysis(x, method = "skip")

Internal missing-values handler. It mirrors ``__missing_values_analysis`` from
pyMannKendall.

When ``method == "skip"`` every observation that is ``NaN`` is removed:

- for a vector, all ``NaN`` entries are dropped;
- for a matrix, every row that contains at least one ``NaN`` is dropped.

Returns ``(x, n)`` where ``n`` is the number of remaining observations
(vector length, or number of matrix rows).
"""
function __missing_values_analysis(x, method = "skip")
    if lowercase(method) == "skip"
        if x isa AbstractVector
            x = x[.!isnan.(x)]
        elseif x isa AbstractMatrix
            keep = Bool[]
            for r in 1:size(x, 1)
                row_has_nan = false
                for c in 1:size(x, 2)
                    if isnan(x[r, c])
                        row_has_nan = true
                        break
                    end
                end
                push!(keep, .!row_has_nan)
            end
            x = x[keep, :]
        end
    end

    n = x isa AbstractVector ? length(x) : size(x, 1)
    return (x, n)
end

"""
    __acf(x, nlags)

Autocorrelation function using the *biased* estimator. It mirrors ``__acf`` from
pyMannKendall.

For lags ``l = 0, 1, ..., nlags`` the autocovariance is

    acov(l) = (1 / n) * sum_{t} (x[t] - mean(x)) * (x[t + l] - mean(x))

and the returned autocorrelation is ``acov(l) / acov(0)``.

The autocorrelation at lag ``l`` is stored at index ``l + 1`` (Julia uses
1-based indexing).
"""
function __acf(x, nlags)
    n = length(x)
    mean_x = sum(x) / n

    # demeaned series
    y = Vector{Float64}(undef, n)
    for i in 1:n
        y[i] = x[i] - mean_x
    end

    # biased autocovariance for lags 0 .. n-1
    acov = zeros(n)
    for l in 0:(n - 1)
        s = 0.0
        for t in 1:(n - l)
            s += y[t] * y[t + l]
        end
        acov[l + 1] = s / n
    end

    m = min(nlags + 1, n)
    if acov[1] != 0.0
        return acov[1:m] ./ acov[1]
    else
        return acov[1:m]
    end
end

"""
    __mk_score(x, n)

Mann-Kendall score ``S``. It mirrors ``__mk_score`` from pyMannKendall.

    S = sum_{i<j} sign(x[j] - x[i])

computed with explicit nested loops.
"""
function __mk_score(x, n)
    s = 0.0
    for k in 1:(n - 1)
        for j in (k + 1):n
            if x[j] > x[k]
                s += 1.0
            elseif x[j] < x[k]
                s -= 1.0
            end
        end
    end
    return s
end

"""
    __variance_s(x, n)

Variance of the Mann-Kendall score ``S``. It mirrors ``__variance_s`` from
pyMannKendall, handling tied values correctly.

If there are no ties (``n == g`` where ``g`` is the number of unique values):

    var(S) = n*(n-1)*(2n+5) / 18

otherwise each tie group of size ``tp`` contributes ``tp*(tp-1)*(2tp+5)`` which
is subtracted from the numerator.
"""
function __variance_s(x, n)
    unique_x = unique(x)
    g = length(unique_x)

    if n == g
        # no ties
        var_s = (n * (n - 1) * (2 * n + 5)) / 18
    else
        # ties are present
        tp = zeros(g)
        for i in 1:g
            cnt = 0
            for j in 1:n
                if x[j] == unique_x[i]
                    cnt += 1
                end
            end
            tp[i] = Float64(cnt)
        end

        sum_term = 0.0
        for i in 1:g
            sum_term += tp[i] * (tp[i] - 1) * (2 * tp[i] + 5)
        end

        var_s = (n * (n - 1) * (2 * n + 5) - sum_term) / 18
    end

    return var_s
end

"""
    __z_score(s, var_s)

Standardized test statistic ``Z`` with continuity correction. It mirrors
``__z_score`` from pyMannKendall.
"""
function __z_score(s, var_s)
    if s == 0
        z = 0.0
    elseif var_s <= 0.0
        # Match Python's np.sqrt(negative) → NaN behavior
        z = NaN
    elseif s > 0
        z = (s - 1) / sqrt(var_s)
    else
        z = (s + 1) / sqrt(var_s)
    end
    return z
end

"""
    __p_value(z, alpha)

Two-tailed p-value, significance decision and trend direction.

    p = 2 * (1 - Phi(|z|))
    h = |z| > Phi^{-1}(1 - alpha/2)

trend is ``"increasing"``, ``"decreasing"`` or ``"no trend"``. Mirrors
``__p_value`` from pyMannKendall.
"""
function __p_value(z, alpha)
    p = 2 * (1 - normcdf(abs(z)))
    h = abs(z) > norminvcdf(1 - alpha / 2)

    if (z < 0) && h
        trend = "decreasing"
    elseif (z > 0) && h
        trend = "increasing"
    else
        trend = "no trend"
    end

    return (p = p, h = h, trend = trend)
end

"""
    __R(x)

Ranks of ``x`` with average ranks for ties. It mirrors ``__R`` from
pyMannKendall (equivalent to ``scipy.stats.rankdata`` with the default
'average' method).

For each element ``x[j]``:

    R[j] = (n + 1 + sum_i sign(x[j] - x[i])) / 2
"""
function __R(x)
    n = length(x)
    R = zeros(n)
    for j in 1:n
        s = 0.0
        for i in 1:n
            if x[j] > x[i]
                s += 1.0
            elseif x[j] < x[i]
                s -= 1.0
            end
        end
        R[j] = (n + 1 + s) / 2
    end
    return R
end

"""
    __K(x, z)

Kendall concordance between ``x`` and ``z``. It mirrors ``__K`` from
pyMannKendall.

    K = sum_{i<j} sign((x[j] - x[i]) * (z[j] - z[i]))
"""
function __K(x, z)
    n = length(x)
    K = 0.0
    for i in 1:(n - 1)
        for j in (i + 1):n
            prod = (x[j] - x[i]) * (z[j] - z[i])
            if prod > 0
                K += 1.0
            elseif prod < 0
                K -= 1.0
            end
        end
    end
    return K
end

"""
    __sens_estimator(x)

All pairwise Theil-Sen slopes ``(x[j] - x[i]) / (j - i)`` for ``i < j``. It
mirrors ``__sens_estimator`` from pyMannKendall.
"""
function __sens_estimator(x)
    n = length(x)
    m = n * (n - 1) ÷ 2
    d = Vector{Float64}(undef, m)
    idx = 1
    for i in 1:(n - 1)
        for j in (i + 1):n
            d[idx] = (x[j] - x[i]) / (j - i)
            idx += 1
        end
    end
    return d
end

"""
    nanmedian(v)

Median of ``v`` ignoring ``NaN`` entries (equivalent to ``numpy.nanmedian``).
"""
function nanmedian(v)
    w = Float64[]
    for val in v
        if !isnan(val)
            push!(w, val)
        end
    end
    return median(w)
end

"""
    __median_valid_index(x)

Median of the 0-based indices of the non-``NaN`` entries of ``x``. This is what
``numpy.median(numpy.arange(n)[~numpy.isnan(x)])`` computes, and it is used for
the intercept of the Kendall-Theil robust line.
"""
function __median_valid_index(x)
    idx = Int[]
    for k in 1:length(x)
        if !isnan(x[k])
            push!(idx, k - 1)
        end
    end
    return median(idx)
end
