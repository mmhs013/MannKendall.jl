# Examples — MannKendall.jl

These examples show how to use **MannKendall.jl** on its own and together with the Julia data-science ecosystem (`DataFrames`, `CSV`, `Statistics`, `Plots`).

All results match [pyMannKendall](https://github.com/mmhs013/pyMannKendall) exactly.

## Setup

```julia
using Pkg
Pkg.add(["MannKendall", "DataFrames", "CSV", "Plots", "Statistics"])

using MannKendall
using DataFrames
using CSV
using Statistics
using Plots
```

---

## 1. Original Mann-Kendall test (vector input)

```julia
# Increasing trend: 1:100
x = collect(1.0:100.0)

result = mk_original_test(x)
println(result)
# (trend = "increasing", h = true, p = 0.0, z = 14.73869998208331,
#  τ = 1.0, s = 4950.0, var_s = 112750.0, slope = 1.0, intercept = 1.0)

# Access fields by name
println("Trend: ", result.trend)   # increasing
println("p-value: ", result.p)     # 0.0
println("Slope: ", result.slope)   # 1.0

# Or unpack directly
trend, h, p, z, τ, s, var_s, slope, intercept = mk_original_test(x)
```

---

## 2. Decreasing and no-trend data

```julia
# Decreasing trend: 100:1
x_dec = collect(100.0:-1.0:1.0)
r_dec = mk_original_test(x_dec)
println(r_dec.trend)  # decreasing
println(r_dec.h)      # true

# No trend: constant data
x_const = ones(100)
r_const = mk_original_test(x_const)
println(r_const.trend)  # no trend
println(r_const.h)      # false
println(r_const.p)      # 1.0
```

---

## 3. Data with ties

```julia
x_ties = [1.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0, 5.0, 5.0, 6.0]
r = mk_original_test(x_ties)
println(r)
# (trend = "increasing", h = true, p = 0.0003568..., z = 3.5701..., τ = 0.8889,
#  s = 40.0, var_s = 119.333..., slope = 0.5, intercept = 0.75)
```

---

## 4. Reading data from CSV with DataFrames

Suppose you have a CSV file `rainfall.csv` with one column of yearly rainfall:

```
year,rainfall
2000,1200.5
2001,1150.2
2002,1230.0
...
```

```julia
df = CSV.read("rainfall.csv", DataFrame)
rainfall = df.rainfall

r = mk_original_test(rainfall)
println("Rainfall trend: ", r.trend)
println("Significant: ", r.h)
println("p-value: ", r.p)
println("Sen's slope: ", r.slope)
```

---

## 5. Multiple stations in one DataFrame (multivariate / regional test)

A DataFrame where each column is a station:

```
station_a,station_b,station_c
1.2,0.8,1.5
1.4,0.9,1.6
...
```

```julia
df = CSV.read("stations.csv", DataFrame)
X = Matrix(df)   # 1000 × 5 matrix, one column per station

# Regional Mann-Kendall test (overall trend across stations)
r_reg = mk_regional_test(X)
println("Regional trend: ", r_reg.trend)
println("Regional p-value: ", r_reg.p)
println("Regional Sen's slope: ", r_reg.slope)

# Multivariate Mann-Kendall test
r_mv = mk_multivariate_test(X)
println("Multivariate trend: ", r_mv.trend)
```

---

## 6. Correlated multivariate test

When the series are correlated, use `mk_correlated_multivariate_test`:

```julia
X = Matrix(df)   # 1000 × 2 or more columns
r_cm = mk_correlated_multivariate_test(X)
println("Correlated multivariate trend: ", r_cm.trend)
println("p-value: ", r_cm.p)
println("Var_s: ", r_cm.var_s)
```

---

## 7. Seasonal Mann-Kendall test (monthly data, period = 12)

```julia
# 120 months = 10 years of monthly data
monthly = rand(120)   # replace with your real monthly series

r_seas = mk_seasonal_test(monthly; period=12)
println("Seasonal trend: ", r_seas.trend)
println("Seasonal p-value: ", r_seas.p)
println("Seasonal Sen's slope: ", r_seas.slope)
println("Seasonal intercept: ", r_seas.intercept)

# Correlated seasonal test (accounts for correlation between months)
r_cs = mk_correlated_seasonal_test(monthly; period=12)
println("Correlated seasonal trend: ", r_cs.trend)
```

---

## 8. Modified tests for serially correlated data

If your data has autocorrelation, use one of the modified tests:

```julia
x = rand(1000)

# Hamed and Rao (1998) variance correction
r_hr = mk_hamed_rao_modification_test(x)
println("Hamed-Rao: ", r_hr.trend, "  p=", r_hr.p)

# Hamed-Rao with only the first 3 significant lags
r_hr3 = mk_hamed_rao_modification_test(x; lag=3)

# Yue and Wang (2004) effective sample size
r_yw = mk_yue_wang_modification_test(x)
r_yw1 = mk_yue_wang_modification_test(x; lag=1)

# Pre-whitening (Yue and Wang 2002)
r_pw = mk_pre_whitening_modification_test(x)

# Trend-free pre-whitening (Yue and Wang 2002)
r_tfpw = mk_trend_free_pre_whitening_modification_test(x)
```

---

## 9. Partial Mann-Kendall test (two columns)

Removes the effect of a covariate from the trend of the response:

```julia
# DataFrame with two columns: response and covariate
df = CSV.read("response_covariate.csv", DataFrame)
X = Matrix(df)   # must have exactly 2 columns

r_partial = mk_partial_test(X)
println("Partial trend: ", r_partial.trend)
println("Partial p-value: ", r_partial.p)
```

---

## 10. Theil-Sen's Slope Estimator

```julia
x = rand(1000)

slope, intercept = sens_slope(x)
println("Theil-Sen slope: ", slope)
println("Intercept: ", intercept)

# Seasonal Theil-Sen (monthly data, period = 12)
slope_s, intercept_s = seasonal_sens_slope(x; period=12)
println("Seasonal slope: ", slope_s)
println("Seasonal intercept: ", intercept_s)
```

---

## 11. Plotting a trend line with Plots.jl

```julia
x = collect(1.0:100.0)
y = 2.0 .* x .+ 5.0 .+ randn(100)   # increasing trend + noise

r = mk_original_test(y)

scatter(x, y, label="data", legend=:topleft)
plot!(x, r.intercept .+ r.slope .* x, label="Theil-Sen trend",
      linecolor=:red, linewidth=2, title="Mann-Kendall: $(r.trend) (p=$(round(r.p, digits=4)))")
```

---

## 12. Complete workflow with a DataFrame

```julia
df = CSV.read("climate.csv", DataFrame)

# Add trend results as new columns
results = [mk_original_test(df[:, col]) for col in names(df)]

trends = [r.trend for r in results]
pvals  = [r.p for r in results]
slopes = [r.slope for r in results]

summary_df = DataFrame(
    variable = names(df),
    trend    = trends,
    p_value  = pvals,
    slope    = slopes,
    significant = [r.h for r in results]
)

println(summary_df)

# Save to CSV
CSV.write("trend_summary.csv", summary_df)
```

---

## Link

These examples are linked from the [README](../README.md). For the full function reference and installation instructions, see the README.
