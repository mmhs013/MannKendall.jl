"""
    sens_slope(x)

Theil-Sen slope estimator (Theil 1950, Sen 1968). Estimates the magnitude of a
monotonic trend. The intercept is computed with the Conover (1980) method.

This mirrors ``sens_slope`` from pyMannKendall.

# Input
- ``x``: a one-dimensional data vector.

# Output
A named tuple ``(slope, intercept)``.

# Example
```julia
using MannKendall
x = rand(120)
slope, intercept = sens_slope(x)
```
"""
function sens_slope(x)
    x, c = __preprocessing(x)
    n = length(x)
    slope = nanmedian(__sens_estimator(x))
    intercept = nanmedian(x) - __median_valid_index(x) * slope
    return (slope = slope, intercept = intercept)
end

"""
    seasonal_sens_slope(x_old; period = 12)

Seasonal Theil-Sen slope estimator (Hipel 1994). Estimates the magnitude of a
monotonic trend when data have seasonal effects. The intercept is computed with
the Conover (1980) method; a full period cycle is considered as one time step.

This mirrors ``seasonal_sens_slope`` from pyMannKendall.

# Input
- ``x_old``: a data vector.
- ``period``: seasonal cycle (12 for monthly, 52 for weekly; default 12).

# Output
A named tuple ``(slope, intercept)``.
"""
function seasonal_sens_slope(x_old; period = 12)
    x, c = __preprocessing(x_old)
    n = length(x)

    if x isa AbstractVector
        if mod(n, period) != 0
            # pad with NaN so that the length is a multiple of period
            pad = period - mod(n, period)
            x = vcat(x, fill(NaN, pad))
        end
        n = length(x)
        # Build the seasonal matrix in numpy (row-major) reshape order:
        # column s holds every period-th element x[s], x[s+period], x[s+2*period], ...
        # (Julia reshape is column-major, so we cannot use it directly here.)
        nrows = Int(n ÷ period)
        xmat = Matrix{Float64}(undef, nrows, period)
        for row in 1:nrows
            for s in 1:period
                xmat[row, s] = x[s + (row - 1) * period]
            end
        end
    else
        xmat = x
    end

    d = Float64[]
    for i in 1:period
        col = xmat[:, i]
        append!(d, __sens_estimator(col))
    end

    slope = nanmedian(d)
    intercept = nanmedian(x_old) - __median_valid_index(x_old) / period * slope

    return (slope = slope, intercept = intercept)
end

"""
    mk_original_test(x; α = 0.05)

Original Mann-Kendall test (Mann 1945, Kendall 1975, Gilbert 1987).

This mirrors ``original_test`` from pyMannKendall.

# Input
- ``x``: a data vector.
- ``α``: significance level (default 0.05).

# Output
A named tuple with fields:
``(trend, h, p, z, τ, s, var_s, slope, intercept)``.
"""
function mk_original_test(x; α = 0.05)
    x, c = __preprocessing(x)
    x, n = __missing_values_analysis(x, "skip")

    s = __mk_score(x, n)
    var_s = __variance_s(x, n)
    τ = s / (0.5 * n * (n - 1))

    z = __z_score(s, var_s)
    pv = __p_value(z, α)
    slope, intercept = sens_slope(x)

    return (
        trend = pv.trend,
        h = pv.h,
        p = pv.p,
        z = z,
        τ = τ,
        s = s,
        var_s = var_s,
        slope = slope,
        intercept = intercept,
    )
end

"""
    mk_hamed_rao_modification_test(x; α = 0.05, lag = nothing)

Modified Mann-Kendall test using the Hamed and Rao (1998) variance correction
to address serial autocorrelation. A variance correction factor ``n_ns``
inflates ``var_s`` based on the significant autocorrelation lags.

This mirrors ``hamed_rao_modification_test`` from pyMannKendall.

# Input
- ``x``: a data vector.
- ``α``: significance level (default 0.05).
- ``lag``: number of first significant lags (default ``nothing``, i.e. all
  lags up to ``n``). Pass e.g. ``lag = 3`` to consider the first 3 lags.

# Output
A named tuple ``(trend, h, p, z, τ, s, var_s, slope, intercept)``.
"""
function mk_hamed_rao_modification_test(x; α = 0.05, lag = nothing)
    x, c = __preprocessing(x)
    x, n = __missing_values_analysis(x, "skip")

    s = __mk_score(x, n)
    var_s = __variance_s(x, n)
    τ = s / (0.5 * n * (n - 1))

    # Hamed and Rao (1998) variance correction
    if lag === nothing
        lag_val = n
    else
        lag_val = lag + 1
    end

    # detrending
    slope, intercept = sens_slope(x)
    x_detrend = x .- collect(1:n) .* slope
    I = __R(x_detrend)

    # account for autocorrelation
    acf_1 = __acf(I, lag_val - 1)
    interval = norminvcdf(1 - α / 2) / sqrt(n)
    upper_bound = 0.0 + interval
    lower_bound = 0.0 - interval

    sni = 0.0
    for i in 1:(lag_val - 1)
        if (acf_1[i + 1] <= upper_bound && acf_1[i + 1] >= lower_bound)
            sni = sni
        else
            sni += (n - i) * (n - i - 1) * (n - i - 2) * acf_1[i + 1]
        end
    end

    n_ns = 1 + (2 / (n * (n - 1) * (n - 2))) * sni
    var_s = var_s * n_ns

    z = __z_score(s, var_s)
    pv = __p_value(z, α)

    return (
        trend = pv.trend,
        h = pv.h,
        p = pv.p,
        z = z,
        τ = τ,
        s = s,
        var_s = var_s,
        slope = slope,
        intercept = intercept,
    )
end

"""
    mk_yue_wang_modification_test(x; α = 0.05, lag = nothing)

Modified Mann-Kendall test using the Yue and Wang (2004) variance correction
(effective sample size) for serially correlated data.

This mirrors ``yue_wang_modification_test`` from pyMannKendall.

# Input
- ``x``: a data vector.
- ``α``: significance level (default 0.05).
- ``lag``: number of first significant lags (default ``nothing``, i.e. all
  lags up to ``n``). Pass e.g. ``lag = 1`` to consider the first 1 lag.

# Output
A named tuple ``(trend, h, p, z, τ, s, var_s, slope, intercept)``.
"""
function mk_yue_wang_modification_test(x; α = 0.05, lag = nothing)
    x, c = __preprocessing(x)
    x, n = __missing_values_analysis(x, "skip")

    s = __mk_score(x, n)
    var_s = __variance_s(x, n)
    τ = s / (0.5 * n * (n - 1))

    # Yue and Wang (2004) variance correction
    if lag === nothing
        lag_val = n
    else
        lag_val = lag + 1
    end

    # detrending
    slope, intercept = sens_slope(x)
    x_detrend = x .- collect(1:n) .* slope

    # account for autocorrelation
    acf_1 = __acf(x_detrend, lag_val - 1)
    idx = collect(1:(lag_val - 1))
    sni = 0.0
    for i in 1:(lag_val - 1)
        sni += (1 - idx[i] / n) * acf_1[i + 1]
    end

    n_ns = 1 + 2 * sni
    var_s = var_s * n_ns

    z = __z_score(s, var_s)
    pv = __p_value(z, α)

    return (
        trend = pv.trend,
        h = pv.h,
        p = pv.p,
        z = z,
        τ = τ,
        s = s,
        var_s = var_s,
        slope = slope,
        intercept = intercept,
    )
end

"""
    mk_pre_whitening_modification_test(x; α = 0.05)

Modified Mann-Kendall test using the Pre-Whitening method (Yue and Wang 2002).
The series is pre-whitened by removing the lag-1 autocorrelation before the
trend test is applied.

This mirrors ``pre_whitening_modification_test`` from pyMannKendall.

# Input
- ``x``: a data vector.
- ``α``: significance level (default 0.05).

# Output
A named tuple ``(trend, h, p, z, τ, s, var_s, slope, intercept)``.
"""
function mk_pre_whitening_modification_test(x; α = 0.05)
    x, c = __preprocessing(x)
    x, n = __missing_values_analysis(x, "skip")

    # PreWhitening
    acf_1 = __acf(x, 1)[2]
    a = collect(0:(n - 2))
    b = collect(1:(n - 1))
    x_pw = Vector{Float64}(undef, n - 1)
    for k in 1:(n - 1)
        x_pw[k] = x[b[k] + 1] - x[a[k] + 1] * acf_1
    end
    n = length(x_pw)

    s = __mk_score(x_pw, n)
    var_s = __variance_s(x_pw, n)
    τ = s / (0.5 * n * (n - 1))

    z = __z_score(s, var_s)
    pv = __p_value(z, α)
    slope, intercept = sens_slope(x)

    return (
        trend = pv.trend,
        h = pv.h,
        p = pv.p,
        z = z,
        τ = τ,
        s = s,
        var_s = var_s,
        slope = slope,
        intercept = intercept,
    )
end

"""
    mk_trend_free_pre_whitening_modification_test(x; α = 0.05)

Modified Mann-Kendall test using the trend-free Pre-Whitening method (Yue and
Wang 2002). The trend component is removed first, then the series is
pre-whitened, and finally the trend is added back before applying the MK test.

This mirrors ``trend_free_pre_whitening_modification_test`` from pyMannKendall.

# Input
- ``x``: a data vector.
- ``α``: significance level (default 0.05).

# Output
A named tuple ``(trend, h, p, z, τ, s, var_s, slope, intercept)``.
"""
function mk_trend_free_pre_whitening_modification_test(x; α = 0.05)
    x, c = __preprocessing(x)
    x, n = __missing_values_analysis(x, "skip")

    # detrending
    slope, intercept = sens_slope(x)
    x_detrend = x .- collect(1:n) .* slope

    # PreWhitening
    acf_1 = __acf(x_detrend, 1)[2]
    a = collect(0:(n - 2))
    b = collect(1:(n - 1))
    x_pw = Vector{Float64}(undef, n - 1)
    for k in 1:(n - 1)
        x_pw[k] = x_detrend[b[k] + 1] - x_detrend[a[k] + 1] * acf_1
    end
    n = length(x_pw)
    x_pw = x_pw .+ collect(1:n) .* slope

    s = __mk_score(x_pw, n)
    var_s = __variance_s(x_pw, n)
    τ = s / (0.5 * n * (n - 1))

    z = __z_score(s, var_s)
    pv = __p_value(z, α)
    slope, intercept = sens_slope(x)

    return (
        trend = pv.trend,
        h = pv.h,
        p = pv.p,
        z = z,
        τ = τ,
        s = s,
        var_s = var_s,
        slope = slope,
        intercept = intercept,
    )
end

"""
    mk_multivariate_test(x; α = 0.05)

Multivariate Mann-Kendall test (Hirsch and Slack 1984; Helsel 2006). Sums the
MK scores and variances over all columns (series) and computes an overall
Kendall Tau.

This mirrors ``multivariate_test`` from pyMannKendall.

# Input
- ``x``: a data matrix with one column per series.
- ``α``: significance level (default 0.05).

# Output
A named tuple ``(trend, h, p, z, τ, s, var_s, slope, intercept)``.
"""
function mk_multivariate_test(x; α = 0.05)
    s = 0.0
    var_s = 0.0
    denom = 0.0

    x_mat, c = __preprocessing(x)

    if x_mat isa AbstractVector
        x_mat = reshape(x_mat, :, 1)
        c = 1
    end

    n = 0
    for i in 1:c
        if c == 1
            x_new, n_i = __missing_values_analysis(x_mat, "skip")
        else
            x_new, n_i = __missing_values_analysis(x_mat[:, i], "skip")
        end
        n = n_i
        s += __mk_score(x_new, n)
        var_s += __variance_s(x_new, n)
        denom += (0.5 * n * (n - 1))
    end

    τ = s / denom

    z = __z_score(s, var_s)
    pv = __p_value(z, α)

    slope, intercept = seasonal_sens_slope(x; period = c)

    return (
        trend = pv.trend,
        h = pv.h,
        p = pv.p,
        z = z,
        τ = τ,
        s = s,
        var_s = var_s,
        slope = slope,
        intercept = intercept,
    )
end

"""
    mk_seasonal_test(x_old; period = 12, α = 0.05)

Seasonal Mann-Kendall test (Hirsch, Slack & Smith 1982). Removes the effect of
seasonality by treating each season as a separate series (via the multivariate
test).

This mirrors ``seasonal_test`` from pyMannKendall.

# Input
- ``x_old``: a data vector.
- ``period``: seasonal cycle (12 for monthly, 52 for weekly; default 12).
- ``α``: significance level (default 0.05).

# Output
A named tuple ``(trend, h, p, z, τ, s, var_s, slope, intercept)``.
"""
function mk_seasonal_test(x_old; period = 12, α = 0.05)
    x, c = __preprocessing(x_old)
    n = length(x)

    if x isa AbstractVector
        if mod(n, period) != 0
            pad = period - mod(n, period)
            x = vcat(x, fill(NaN, pad))
        end
        n = length(x)
        # Build the seasonal matrix in numpy (row-major) reshape order:
        # column s holds every period-th element x[s], x[s+period], x[s+2*period], ...
        # (Julia reshape is column-major, so we cannot use it directly here.)
        nrows = Int(n ÷ period)
        x_mat = Matrix{Float64}(undef, nrows, period)
        for row in 1:nrows
            for s in 1:period
                x_mat[row, s] = x[s + (row - 1) * period]
            end
        end
    else
        x_mat = x
    end

    trend, h, p, z, τ, s, var_s, slope, intercept =
        mk_multivariate_test(x_mat; α = α)

    return (
        trend = trend,
        h = h,
        p = p,
        z = z,
        τ = τ,
        s = s,
        var_s = var_s,
        slope = slope,
        intercept = intercept,
    )
end

"""
    mk_regional_test(x; α = 0.05)

Regional Mann-Kendall test (Helsel 2006). Computes the overall trend across a
set of stations/locations, where each column of ``x`` is one location.

This mirrors ``regional_test`` from pyMannKendall.

# Input
- ``x``: a data matrix, one column per station.
- ``α``: significance level (default 0.05).

# Output
A named tuple ``(trend, h, p, z, τ, s, var_s, slope, intercept)``.
"""
function mk_regional_test(x; α = 0.05)
    return mk_multivariate_test(x; α = α)
end

"""
    mk_correlated_multivariate_test(x; α = 0.05)

Correlated Multivariate Mann-Kendall test (Libiseller and Grimvall 2002).
Accounts for correlation between the different series/columns of ``x`` when
computing the variance.

This mirrors ``correlated_multivariate_test`` from pyMannKendall.

# Input
- ``x``: a data matrix, one column per series.
- ``α``: significance level (default 0.05).

# Output
A named tuple ``(trend, h, p, z, τ, s, var_s, slope, intercept)``.
"""
function mk_correlated_multivariate_test(x; α = 0.05)
    x_mat, c = __preprocessing(x)
    x_mat, n = __missing_values_analysis(x_mat, "skip")

    s = 0.0
    denom = 0.0

    for i in 1:c
        s += __mk_score(x_mat[:, i], n)
        denom += (0.5 * n * (n - 1))
    end

    τ = s / denom

    Gamma = zeros(c, c)

    for i in 1:c
        for j in 1:c
            k = __K(x_mat[:, i], x_mat[:, j])
            ri = __R(x_mat[:, i])
            rj = __R(x_mat[:, j])
            Gamma[i, j] = (k + 4 * sum(ri .* rj) - n * (n + 1)^2) / 3
        end
    end

    var_s = sum(Gamma)

    z = s / sqrt(var_s)
    pv = __p_value(z, α)
    slope, intercept = seasonal_sens_slope(x_mat; period = c)

    return (
        trend = pv.trend,
        h = pv.h,
        p = pv.p,
        z = z,
        τ = τ,
        s = s,
        var_s = var_s,
        slope = slope,
        intercept = intercept,
    )
end

"""
    mk_correlated_seasonal_test(x_old; period = 12, α = 0.05)

Correlated Seasonal Mann-Kendall test (Hipel 1994). Used when the time series is
significantly correlated with preceding months/seasons.

This mirrors ``correlated_seasonal_test`` from pyMannKendall.

# Input
- ``x_old``: a data vector.
- ``period``: seasonal cycle (default 12).
- ``α``: significance level (default 0.05).

# Output
A named tuple ``(trend, h, p, z, τ, s, var_s, slope, intercept)``.
"""
function mk_correlated_seasonal_test(x_old; period = 12, α = 0.05)
    x, c = __preprocessing(x_old)
    n = length(x)

    if x isa AbstractVector
        if mod(n, period) != 0
            pad = period - mod(n, period)
            x = vcat(x, fill(NaN, pad))
        end
        n = length(x)
        # Build the seasonal matrix in numpy (row-major) reshape order:
        # column s holds every period-th element x[s], x[s+period], x[s+2*period], ...
        # (Julia reshape is column-major, so we cannot use it directly here.)
        nrows = Int(n ÷ period)
        x_mat = Matrix{Float64}(undef, nrows, period)
        for row in 1:nrows
            for s in 1:period
                x_mat[row, s] = x[s + (row - 1) * period]
            end
        end
    else
        x_mat = x
    end

    return mk_correlated_multivariate_test(x_mat; α = α)
end

"""
    mk_partial_test(x; α = 0.05)

Partial Mann-Kendall test (Libiseller and Grimvall 2002). Requires two
parameters (two columns): a response variable and an independent covariate. It
removes the effect of the covariate from the trend estimate of the response.

This mirrors ``partial_test`` from pyMannKendall.

# Input
- ``x``: a data matrix with exactly 2 columns.
- ``α``: significance level (default 0.05).

# Output
A named tuple ``(trend, h, p, z, τ, s, var_s, slope, intercept)``.
"""
function mk_partial_test(x; α = 0.05)
    x_proc, c = __preprocessing(x)
    x_proc, n = __missing_values_analysis(x_proc, "skip")

    if c != 2
        error(
            "Partial Mann Kendall test requires two parameters/columns. " *
            "Here column number $c is not equal to 2.",
        )
    end

    x = x_proc[:, 1]
    y = x_proc[:, 2]

    x_score = __mk_score(x, n)
    y_score = __mk_score(y, n)

    k = __K(x, y)
    rx = __R(x)
    ry = __R(y)

    sigma = (k + 4 * sum(rx .* ry) - n * (n + 1)^2) / 3
    rho = sigma / (n * (n - 1) * (2 * n + 5) / 18)

    s = x_score - rho * y_score
    var_s = (1 - rho^2) * (n * (n - 1) * (2 * n + 5)) / 18

    τ = x_score / (0.5 * n * (n - 1))

    z = s / sqrt(var_s)
    pv = __p_value(z, α)
    slope, intercept = sens_slope(x_proc[:, 1])

    return (
        trend = pv.trend,
        h = pv.h,
        p = pv.p,
        z = z,
        τ = τ,
        s = s,
        var_s = var_s,
        slope = slope,
        intercept = intercept,
    )
end
