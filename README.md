# MannKendall.jl

[![Build Status](https://travis-ci.com/mmhs013/MannKendall.jl.svg?branch=master)](https://travis-ci.com/mmhs013/MannKendall.jl)
[![Julia Registry](https://img.shields.io/badge/Julia-Registered-blue)](https://github.com/JuliaRegistries/General)

`MannKendall.jl` is a pure Julia package for Mann-Kendall trend analysis. It is the Julia version of [pyMannKendall](https://github.com/mmhs013/pyMannKendall), covering **all** functions of the pyMannKendall package.

The Mann-Kendall Trend Test (MK test) is used to analyze time series data for consistently increasing or decreasing trends (monotonic trends). It is a non-parametric test, which means it works for all distributions, but data should have no serial correlation. When serial correlation is present, several modified Mann-Kendall tests are available (Hamed and Rao, Yue and Wang, Pre-Whitening, Trend-free Pre-Whitening, etc.). Seasonal, multivariate, regional, correlated, and partial variants are also included.

## What is implemented?

`MannKendall.jl` brings together almost all types of Mann-Kendall tests. Currently the package has **11 Mann-Kendall tests** and **2 Sen's slope estimator functions**:

1. **Original Mann-Kendall test** (`mk_original_test`) — nonparametric test, no serial correlation or seasonal effects.
2. **Hamed and Rao Modified MK Test** (`mk_hamed_rao_modification_test`) — variance correction for serial autocorrelation (Hamed & Rao 1998).
3. **Yue and Wang Modified MK Test** (`mk_yue_wang_modification_test`) — variance correction / effective sample size (Yue & Wang 2004).
4. **Pre-Whitening Modified MK Test** (`mk_pre_whitening_modification_test`) — pre-whiten the series before testing (Yue & Wang 2002).
5. **Trend-free Pre-Whitening Modified MK Test** (`mk_trend_free_pre_whitening_modification_test`) — remove trend, pre-whiten, add trend back (Yue & Wang 2002).
6. **Multivariate MK Test** (`mk_multivariate_test`) — MK test for multiple parameters (Hirsch & Slack 1984).
7. **Seasonal MK Test** (`mk_seasonal_test`) — seasonal trend for monthly/weekly data (Hirsch, Slack & Smith 1982).
8. **Regional MK Test** (`mk_regional_test`) — overall trend across stations/locations (Helsel 2006).
9. **Correlated Multivariate MK Test** (`mk_correlated_multivariate_test`) — correlated parameters (Libiseller & Grimvall 2002).
10. **Correlated Seasonal MK Test** (`mk_correlated_seasonal_test`) — seasonal data with correlation between months (Hipel 1994).
11. **Partial MK Test** (`mk_partial_test`) — removes the effect of a covariate; requires two columns (Libiseller & Grimvall 2002).
12. **Theil-Sen's Slope Estimator** (`sens_slope`) — magnitude of monotonic trend (Theil 1950, Sen 1968).
13. **Seasonal Theil-Sen's Slope Estimator** (`seasonal_sens_slope`) — magnitude when data has seasonal effects (Hipel 1994).

> All results have been verified to match [pyMannKendall](https://github.com/mmhs013/pyMannKendall) exactly.

## Function details

All Mann-Kendall test functions share similar input arguments:

- **`x`**: a data vector (or matrix for multivariate/regional/correlated/partial tests)
- **`α`**: significance level (default `0.05`)
- **`lag`**: number of first significant lags (`hamed_rao_modification_test` and `yue_wang_modification_test` only)
- **`period`**: seasonal cycle — 12 for monthly, 52 for weekly (default `12`; seasonal tests only)

All Mann-Kendall tests return a named tuple containing:

- **`trend`**: `"increasing"`, `"decreasing"`, or `"no trend"`
- **`h`**: `true` if a trend is present, `false` if absent
- **`p`**: p-value of the significance test
- **`z`**: normalized test statistic
- **`τ`** (`tau`): Kendall Tau
- **`s`**: Mann-Kendall's score
- **`var_s`**: Variance of S
- **`slope`**: Theil-Sen estimator/slope
- **`intercept`**: intercept of the Kendall-Theil robust line (for seasonal tests, a full period cycle is one time step)

The Sen's slope functions return only `(slope, intercept)`.

## Examples

Several examples, including how to use the package with [DataFrames](https://github.com/JuliaData/DataFrames.jl) and [CSV](https://github.com/JuliaData/CSV.jl), are provided [here](example.md).

## Installation

You can install `MannKendall.jl` from the Julia package manager:

```julia
using Pkg
Pkg.add("MannKendall")
```

Or, to use the latest development version from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/mmhs013/MannKendall.jl")
```

## Usage

A quick example:

```julia
using MannKendall

# Data generation for analysis
x = collect(1.0:100.0)

result = mk_original_test(x)
println(result)
```

Output:

```
(trend = "increasing", h = true, p = 0.0, z = 14.73869998208331, τ = 1.0, s = 4950.0, var_s = 112750.0, slope = 1.0, intercept = 1.0)
```

The output is a named tuple, so you can access fields by name:

```julia
println(result.slope)
```

Or unpack directly:

```julia
trend, h, p, z, τ, s, var_s, slope, intercept = mk_original_test(x)
```

## Running the tests

```julia
using Pkg
Pkg.test("MannKendall")
```

## Citation

If you publish results for which you used `MannKendall.jl` / `pyMannKendall`, please cite [Hussain et al., (2019)](https://doi.org/10.21105/joss.01556):

> Hussain et al., (2019). pyMannKendall: a python package for non parametric Mann Kendall family of trend tests.. Journal of Open Source Software, 4(39), 1556, https://doi.org/10.21105/joss.01556

BibTeX:

```bibtex
@article{Hussain2019pyMannKendall,
	journal = {Journal of Open Source Software},
	doi = {10.21105/joss.01556},
	issn = {2475-9066},
	number = {39},
	publisher = {The Open Journal},
	title = {pyMannKendall: a python package for non parametric Mann Kendall family of trend tests.},
	url = {http://dx.doi.org/10.21105/joss.01556},
	volume = {4},
	author = {Hussain, Md. and Mahmud, Ishtiak},
	pages = {1556},
	date = {2019-07-25},
	year = {2019},
	month = {7},
	day = {25},
}
```

## Contributions

`MannKendall.jl` is a community project and welcomes contributions.

## References

1. Bari, S. H., Rahman, M. T. U., Hoque, M. A., & Hussain, M. M. (2016). Analysis of seasonal and annual rainfall trends in the northern region of Bangladesh. *Atmospheric Research*, 176, 148-158.
2. Conover, W.J. (1980). *Practical Nonparametric Statistics* (2nd Ed.). John Wiley and Sons.
3. Hamed, K. H., & Rao, A. R. (1998). A modified Mann-Kendall trend test for autocorrelated data. *Journal of Hydrology*, 204(1-4), 182-196.
4. Helsel, D. R., & Frans, L. M. (2006). Regional Kendall test for trend. *Environmental Science & Technology*, 40(13), 4066-4073.
5. Hipel, K. W., & McLeod, A. I. (1994). *Time Series Modelling of Water Resources and Environmental Systems* (Vol. 45). Elsevier.
6. Hirsch, R. M., Slack, J. R., & Smith, R. A. (1982). Techniques of trend analysis for monthly water quality data. *Water Resources Research*, 18(1), 107-121.
7. Kendall, M. (1975). *Rank Correlation Measures*. Charles Griffin, London.
8. Libiseller, C., & Grimvall, A. (2002). Performance of partial Mann-Kendall tests for trend detection in the presence of covariates. *Environmetrics*, 13(1), 71-84.
9. Mann, H. B. (1945). Nonparametric tests against trend. *Econometrica*, 245-259.
10. Sen, P. K. (1968). Estimates of the regression coefficient based on Kendall's tau. *Journal of the American Statistical Association*, 63(324), 1379-1389.
11. Theil, H. (1950). A rank-invariant method of linear and polynomial regression analysis. *Indagationes Mathematicae*, 53, 1397-1412.
12. Yue, S., & Wang, C. (2004). The Mann-Kendall test modified by effective sample size to detect trend in serially correlated hydrological series. *Water Resources Management*, 18(3), 201-218.
13. Yue, S., & Wang, C. Y. (2002). Applicability of prewhitening to eliminate the influence of serial correlation on the Mann-Kendall test. *Water Resources Research*, 38(6).
14. Yue, S., Pilon, P., Phinney, B., & Cavadias, G. (2002). The influence of autocorrelation on the ability to detect trend in hydrological series. *Hydrological Processes*, 16(9), 1807-1829.
