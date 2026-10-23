using MannKendall
using Test
using Statistics

# -----------------------------------------------------------------------------
# Test fixtures (mirrors test_pymannkendall.py from pyMannKendall)
# -----------------------------------------------------------------------------

# NoTrendData: 360 random values all the same number
const NoTrendData = ones(360) * Float64(rand(0:9))

# NoTrend2dData: 360×2 matrix, all same number
const NoTrend2dData = ones(360, 2) * Float64(rand(0:9))

# TrendData: 360 trendy data with slope ~ 1
const TrendData = collect(0.0:359.0) .+ rand(360) ./ 10^13

# arbitrary_1d_data (359 values, includes NaNs) — exact data from test_pymannkendall.py
const arbitrary_1d_data = Float64[
    32.,  20.,  25., 189., 240., 193., 379., 278., 301.,   0.,   0.,
    82.,   0.,   4.,  NaN,  NaN, 121., 234., 360., 262., 120.,  30.,
    11.,   1.,   7.,   3.,  31.,  31., 355., 102., 248., 274., 308.,
    NaN,   5.,  26.,  11.,  16.,   6.,  48., 388., 539., 431., 272.,
   404., 186.,   0.,   2.,   0.,   4.,   1.,  54., 272., 459., 235.,
   164., 365., 135.,   2.,  NaN,  NaN,   4.,   0., 128., 210., 163.,
   446., 225., 462., 467.,  19.,  13.,   0.,   3.,  17., 132., 178.,
   338., 525., 623., 145.,  31.,  19.,   3.,   0.,  29.,  25.,  87.,
   259., 756., 486., 180.,  43.,  92.,   1.,   0.,  16.,   2.,   0.,
   130., 253., 594., 111., 273.,  30.,   0.,   4.,   0.,  27.,  24.,
    41., 292., 378., 499., 265., 320., 227.,   4.,   0.,   4.,  14.,
     8.,  48., 416., 240., 404., 207., 733., 105.,   0., 112.,   0.,
    14.,   0.,  30., 140., 202., 289., 159., 424., 106.,   3.,
     0.,  65.,   3.,  14.,  58., 268., 466., 432., 266., 240.,  95.,
     1.,   0.,  10.,  26.,   4., 114.,  94., 289., 173., 208., 263.,
   156.,   5.,   0.,  16.,  16.,  14.,   0., 111., 475., 534., 432.,
   471., 117.,  70.,   1.,   3.,  28.,   7., 401., 184., 283., 338.,
   171., 335., 176.,   0.,   0.,  10.,  11.,   9., 140., 102., 208.,
   298., 245., 220.,  29.,   2.,  27.,  10.,  13.,  26.,  84., 143.,
   367., 749., 563., 283.,  10.,   0.,   0.,   0.,   0.,   9.,
   246., 265., 343., 429., 168., 133.,  17.,   0.,  18.,  35.,  76.,
   158., 272., 250., 190., 289., 466.,  84.,   0.,   0.,   0.,   0.,
     0.,   0.,  22., 217., 299., 185., 115., 344., 203.,   8.,  NaN,  NaN,
     0.,   5., 284., 123., 254., 476., 496., 326.,  27.,  20.,   0.,
     4.,  53.,  72., 113., 214., 364., 219., 220., 156., 264.,   0.,
    13.,   0.,   0.,  45.,  90., 137., 638., 529., 261., 206., 251.,
     0.,   0.,   5.,   9.,  58.,  72., 138., 130., 471., 328., 356.,
   523.,   0.,   1.,   0.,   0.,  12., 143., 193., 184., 192., 138.,
   174.,  69.,   1.,   0.,   0.,  18.,  25.,  28.,  92., 732., 320.,
   256., 302., 131.,  15.,   0.,  27.,   0.,  22.,  20., 213., 393.,
   474., 374., 109., 159.,   0.,   0.,   0.,   3.,   3.,  49., 205.,
   128., 194., 570., 169.,  89.,   0.,   0.,   0.,   0.,   0.,  26.,
   185., 286.,  92., 225., 244., 190.,   3.,  20.]

# arbitrary_2d_data (80×2 matrix, includes NaNs)
const arbitrary_2d_data = Float64[
    490.  458.; 540.  469.; 220. 4630.; 390.  321.; 450.  541.;
   230. 1640.; 360. 1060.; 460.  264.; 430.  665.; 430.  680.;
   620.  650.; 460.  NaN ; 450.  380.; 580.  325.; 350. 1020.;
   440.  460.; 530.  583.; 380.  777.; 440. 1230.; 430.  565.;
   680.  533.; 250.  493.; NaN  3810.; 450.  469.; 500.  473.;
   510.  593.; 490.  500.; 700.  266.; 420.  495.; 710.  245.;
   430.  736.; 410.  508.; 700.  578.; 260.  459.; 260. 4670.;
   500.  503.; 450.  469.; 500.  314.; 620.  432.; 670.  279.;
   NaN   542.; 470.  499.; 370.  741.; 410.  569.; 540.  360.;
   550.  513.; 220. 3910.; 460.  364.; 390.  472.; 550.  245.;
   320.  NaN ; 570.  224.; 480.  342.; 520.  732.; 620.  240.;
   520.  472.; 430.  679.; 400. 1080.; 430.  920.; 490.  488.;
   560.  NaN ; 370.  595.; 460.  295.; 390.  542.; 330. 1500.;
   350. 1080.; 480.  334.; 390.  423.; 500.  216.; 410.  366.;
   470.  750.; 280. 1260.; 510.  223.; NaN   462.; 310. 7640.;
   230. 2340.; 470.  239.; 330. 1400.; 320. 3070.; 500.  244.]

# -----------------------------------------------------------------------------
# sens_slope
# -----------------------------------------------------------------------------
@testset "sens_slope" begin
    # no trend data
    NoTrendRes = sens_slope(NoTrendData)
    @test NoTrendRes.slope == 0.0

    # trendy data
    TrendRes = sens_slope(TrendData)
    @test TrendRes.slope == 1.0
    @test round(TrendRes.intercept) == 0.0

    # arbitrary 1d data
    result = sens_slope(arbitrary_1d_data)
    @test result.slope == -0.00586515307877537
    @test result.intercept == 94.05866013071895
end

# -----------------------------------------------------------------------------
# seasonal_sens_slope
# -----------------------------------------------------------------------------
@testset "seasonal_sens_slope" begin
    # no trend data
    NoTrendRes = seasonal_sens_slope(NoTrendData)
    @test NoTrendRes.slope == 0.0

    # trendy data
    TrendRes = seasonal_sens_slope(TrendData)
    @test TrendRes.slope == 12.0
    @test round(TrendRes.intercept) == 0.0

    # arbitrary 1d data
    result = seasonal_sens_slope(arbitrary_1d_data)
    @test result.slope == 0.0
    @test result.intercept == 93.0
end

# -----------------------------------------------------------------------------
# mk_original_test
# -----------------------------------------------------------------------------
@testset "mk_original_test" begin
    # no trend data
    NoTrendRes = mk_original_test(NoTrendData)
    @test NoTrendRes.trend == "no trend"
    @test NoTrendRes.h == false
    @test NoTrendRes.p == 1.0
    @test NoTrendRes.z == 0
    @test NoTrendRes.τ == 0.0
    @test NoTrendRes.s == 0.0
    @test NoTrendRes.var_s == 0.0

    # trendy data
    TrendRes = mk_original_test(TrendData)
    @test TrendRes.trend == "increasing"
    @test TrendRes.h == true
    @test TrendRes.p == 0.0
    @test TrendRes.τ == 1.0
    @test TrendRes.s == 64620.0

    # arbitrary data
    result = mk_original_test(arbitrary_1d_data)
    @test result.trend == "no trend"
    @test result.h == false
    @test isapprox(result.p, 0.3711454868289752)
    @test isapprox(result.z, -0.8943297502164035)
    @test isapprox(result.τ, -0.03188940688940689)
    @test result.s == -1970.0
    @test isapprox(result.var_s, 4847257.333333333)
end

# -----------------------------------------------------------------------------
# mk_hamed_rao_modification_test
# -----------------------------------------------------------------------------
@testset "mk_hamed_rao_modification_test" begin
    # no trend data
    NoTrendRes = mk_hamed_rao_modification_test(NoTrendData)
    @test NoTrendRes.trend == "no trend"
    @test NoTrendRes.h == false
    @test NoTrendRes.p == 1.0
    @test NoTrendRes.z == 0
    @test NoTrendRes.τ == 0.0
    @test NoTrendRes.s == 0.0

    # trendy data
    TrendRes = mk_hamed_rao_modification_test(TrendData)
    @test TrendRes.trend == "increasing"
    @test TrendRes.h == true
    @test TrendRes.p == 0.0
    @test TrendRes.τ == 1.0
    @test TrendRes.s == 64620.0

    # arbitrary data — pyMannKendall returns NaN p/z and negative var_s for this data
    result = mk_hamed_rao_modification_test(arbitrary_1d_data)
    @test result.trend == "no trend"
    @test isnan(result.p)
    @test isnan(result.z)
    @test result.τ == -0.03188940688940689
    @test result.s == -1970.0
    @test result.var_s < 0  # negative variance due to Hamed-Rao correction
end

# -----------------------------------------------------------------------------
# mk_hamed_rao_modification_test (lag=3)
# -----------------------------------------------------------------------------
@testset "mk_hamed_rao_modification_test_lag3" begin
    # no trend data
    NoTrendRes = mk_hamed_rao_modification_test(NoTrendData; lag=3)
    @test NoTrendRes.trend == "no trend"
    @test NoTrendRes.h == false
    @test NoTrendRes.p == 1.0
    @test NoTrendRes.z == 0
    @test NoTrendRes.τ == 0.0
    @test NoTrendRes.s == 0.0

    # trendy data
    TrendRes = mk_hamed_rao_modification_test(TrendData; lag=3)
    @test TrendRes.trend == "increasing"
    @test TrendRes.h == true
    @test TrendRes.p == 0.0
    @test TrendRes.τ == 1.0
    @test TrendRes.s == 64620.0

    # arbitrary data
    result = mk_hamed_rao_modification_test(arbitrary_1d_data; lag=3)
    @test result.trend == "no trend"
    @test result.h == false
    @test isapprox(result.p, 0.5998889143127422; rtol=1e-3)
    @test isapprox(result.z, -0.5245602663760147; rtol=1e-3)
    @test isapprox(result.τ, -0.03188940688940689)
    @test result.s == -1970.0
    @test isapprox(result.var_s, 14089664.434304914; rtol=1e-3)
end

# -----------------------------------------------------------------------------
# mk_yue_wang_modification_test
# -----------------------------------------------------------------------------
@testset "mk_yue_wang_modification_test" begin
    # no trend data
    NoTrendRes = mk_yue_wang_modification_test(NoTrendData)
    @test NoTrendRes.trend == "no trend"
    @test NoTrendRes.h == false
    @test NoTrendRes.p == 1.0
    @test NoTrendRes.z == 0
    @test NoTrendRes.τ == 0.0
    @test NoTrendRes.s == 0.0

    # trendy data
    TrendRes = mk_yue_wang_modification_test(TrendData)
    @test TrendRes.trend == "increasing"
    @test TrendRes.h == true
    @test TrendRes.p == 0.0
    @test TrendRes.τ == 1.0
    @test TrendRes.s == 64620.0

    # arbitrary data
    result = mk_yue_wang_modification_test(arbitrary_1d_data)
    @test result.trend == "decreasing"
    @test result.h == true
    @test isapprox(result.p, 0.00662187451423879; rtol=1e-2)
    @test isapprox(result.z, -2.715285060390629; rtol=1e-2)
    @test isapprox(result.τ, -0.03188940688940689)
    @test result.s == -1970.0
    @test isapprox(result.var_s, 525848.4184261264; rtol=1e-2)
end

# -----------------------------------------------------------------------------
# mk_yue_wang_modification_test (lag=1)
# -----------------------------------------------------------------------------
@testset "mk_yue_wang_modification_test_lag1" begin
    # no trend data
    NoTrendRes = mk_yue_wang_modification_test(NoTrendData; lag=1)
    @test NoTrendRes.trend == "no trend"
    @test NoTrendRes.h == false
    @test NoTrendRes.p == 1.0
    @test NoTrendRes.z == 0
    @test NoTrendRes.τ == 0.0
    @test NoTrendRes.s == 0.0

    # trendy data
    TrendRes = mk_yue_wang_modification_test(TrendData; lag=1)
    @test TrendRes.trend == "increasing"
    @test TrendRes.h == true
    @test TrendRes.p == 0.0
    @test TrendRes.τ == 1.0
    @test TrendRes.s == 64620.0

    # arbitrary data
    result = mk_yue_wang_modification_test(arbitrary_1d_data; lag=1)
    @test result.trend == "no trend"
    @test result.h == false
    @test isapprox(result.p, 0.5399544327386481; rtol=1e-3)
    @test isapprox(result.z, -0.6128818989871667; rtol=1e-3)
    @test isapprox(result.τ, -0.03188940688940689)
    @test result.s == -1970.0
    @test isapprox(result.var_s, 10321381.774215667; rtol=1e-3)
end

# -----------------------------------------------------------------------------
# mk_pre_whitening_modification_test
# -----------------------------------------------------------------------------
@testset "mk_pre_whitening_modification_test" begin
    # no trend data
    NoTrendRes = mk_pre_whitening_modification_test(NoTrendData)
    @test NoTrendRes.trend == "no trend"
    @test NoTrendRes.h == false
    @test NoTrendRes.p == 1.0
    @test NoTrendRes.z == 0
    @test NoTrendRes.τ == 0.0

    # trendy data
    TrendRes = mk_pre_whitening_modification_test(TrendData)
    @test TrendRes.trend == "increasing"
    @test TrendRes.h == true
    @test TrendRes.p == 0.0

    # arbitrary data
    result = mk_pre_whitening_modification_test(arbitrary_1d_data)
    @test result.trend == "no trend"
    @test result.h == false
    @test isapprox(result.p, 0.9390309309037164)
    @test isapprox(result.z, -0.0764879116548233)
    @test isapprox(result.τ, -0.0027513227513227515)
    @test result.s == -169.0
    @test result.var_s == 4824285.0
end

# -----------------------------------------------------------------------------
# mk_trend_free_pre_whitening_modification_test
# -----------------------------------------------------------------------------
@testset "mk_trend_free_pre_whitening_modification_test" begin
    # no trend data
    NoTrendRes = mk_trend_free_pre_whitening_modification_test(NoTrendData)
    @test NoTrendRes.trend == "no trend"
    @test NoTrendRes.h == false
    @test NoTrendRes.p == 1.0
    @test NoTrendRes.z == 0
    @test NoTrendRes.τ == 0.0

    # trendy data
    TrendRes = mk_trend_free_pre_whitening_modification_test(TrendData)
    @test TrendRes.trend == "increasing"
    @test TrendRes.h == true
    @test TrendRes.p == 0.0
    @test TrendRes.τ == 1.0

    # arbitrary data
    result = mk_trend_free_pre_whitening_modification_test(arbitrary_1d_data)
    @test result.trend == "no trend"
    @test result.h == false
    @test isapprox(result.p, 0.7910489075301421; rtol=1e-3)
    @test isapprox(result.z, -0.2649487957519801; rtol=1e-3)
    @test isapprox(result.τ, -0.009491249491249491; rtol=1e-3)
    @test result.s == -583.0
    @test result.var_s == 4825275.0
end

# -----------------------------------------------------------------------------
# mk_seasonal_test
# -----------------------------------------------------------------------------
@testset "mk_seasonal_test" begin
    # no trend data
    NoTrendRes = mk_seasonal_test(NoTrendData; period=12)
    @test NoTrendRes.trend == "no trend"
    @test NoTrendRes.h == false
    @test NoTrendRes.p == 1.0
    @test NoTrendRes.z == 0
    @test NoTrendRes.τ == 0.0
    @test NoTrendRes.s == 0.0

    # trendy data
    TrendRes = mk_seasonal_test(TrendData; period=12)
    @test TrendRes.trend == "increasing"
    @test TrendRes.h == true
    @test TrendRes.p == 0.0
    @test TrendRes.τ == 1.0
    @test TrendRes.s == 5220.0

    # arbitrary data
    result = mk_seasonal_test(arbitrary_1d_data; period=12)
    @test result.trend == "no trend"
    @test result.h == false
    @test isapprox(result.p, 0.3510032664777658)
    @test isapprox(result.z, -0.9326450474155363)
    @test isapprox(result.τ, -0.03486275295531958)
    @test result.s == -174.0
    @test result.var_s == 34408.0
end

# -----------------------------------------------------------------------------
# mk_regional_test
# -----------------------------------------------------------------------------
@testset "mk_regional_test" begin
    # no trend data
    NoTrendRes = mk_regional_test(NoTrend2dData)
    @test NoTrendRes.trend == "no trend"
    @test NoTrendRes.h == false
    @test NoTrendRes.p == 1.0
    @test NoTrendRes.z == 0
    @test NoTrendRes.τ == 0.0
    @test NoTrendRes.s == 0.0
    @test NoTrendRes.var_s == 0.0
    @test NoTrendRes.slope == 0.0

    # arbitrary data
    result = mk_regional_test(arbitrary_2d_data)
    @test result.trend == "no trend"
    @test result.h == false
    @test isapprox(result.p, 0.30890513282283427; rtol=1e-3)
    @test isapprox(result.z, -1.0175220823429532; rtol=1e-3)
    @test isapprox(result.τ, -0.056049213943950786; rtol=1e-3)
    @test result.s == -328.0
    @test result.var_s == 103278.0
    @test isapprox(result.slope, -0.6138147566718996; rtol=1e-3)
    @test isapprox(result.intercept, 495.399136577708; rtol=1e-3)
end

# -----------------------------------------------------------------------------
# mk_correlated_multivariate_test
# -----------------------------------------------------------------------------
@testset "mk_correlated_multivariate_test" begin
    # no trend data
    NoTrendRes = mk_correlated_multivariate_test(NoTrend2dData)
    @test NoTrendRes.trend == "no trend"
    @test NoTrendRes.h == false
    @test NoTrendRes.τ == 0.0
    @test NoTrendRes.s == 0.0
    @test NoTrendRes.var_s == 0.0
    @test NoTrendRes.slope == 0.0

    # arbitrary data
    result = mk_correlated_multivariate_test(arbitrary_2d_data)
    @test result.trend == "no trend"
    @test result.h == false
    @test isapprox(result.p, 0.1161897839407171; rtol=1e-3)
    @test isapprox(result.z, -1.5709692828547697; rtol=1e-3)
    @test isapprox(result.τ, -0.05349870418363569; rtol=1e-3)
    @test result.s == -289.0
    @test isapprox(result.var_s, 33842.333333333336; rtol=1e-3)
    @test isapprox(result.slope, -0.6417295813315032; rtol=1e-3)
    @test isapprox(result.intercept, 494.5835621139327; rtol=1e-3)
end

# -----------------------------------------------------------------------------
# mk_correlated_seasonal_test
# -----------------------------------------------------------------------------
@testset "mk_correlated_seasonal_test" begin
    # no trend data
    NoTrendRes = mk_correlated_seasonal_test(NoTrendData; period=12)
    @test NoTrendRes.trend == "no trend"
    @test NoTrendRes.h == false
    @test NoTrendRes.τ == 0.0
    @test NoTrendRes.s == 0.0

    # trendy data
    TrendRes = mk_correlated_seasonal_test(TrendData; period=12)
    @test TrendRes.trend == "increasing"
    @test TrendRes.h == true
    @test round(TrendRes.p) == 0
    @test TrendRes.τ == 1.0
    @test TrendRes.s == 5220.0

    # arbitrary data
    result = mk_correlated_seasonal_test(arbitrary_1d_data; period=12)
    @test result.trend == "no trend"
    @test result.h == false
    @test isapprox(result.p, 0.15855883533437787)
    @test isapprox(result.z, -1.4099351403890978)
    @test isapprox(result.τ, -0.05495169082125604)
    @test result.s == -182.0
    @test isapprox(result.var_s, 16662.666666666668)
    @test result.slope == 0.0
    @test result.intercept == 91.0
end

# -----------------------------------------------------------------------------
# mk_partial_test
# -----------------------------------------------------------------------------
@testset "mk_partial_test" begin
    # no trend data
    NoTrendRes = mk_partial_test(NoTrend2dData)
    @test NoTrendRes.trend == "no trend"
    @test NoTrendRes.h == false
    @test NoTrendRes.p == 1.0
    @test NoTrendRes.z == 0
    @test NoTrendRes.τ == 0.0
    @test NoTrendRes.s == 0.0
    @test NoTrendRes.var_s == 5205500.0

    # arbitrary data
    result = mk_partial_test(arbitrary_2d_data)
    @test result.trend == "no trend"
    @test result.h == false
    @test isapprox(result.p, 0.1214998019842275; rtol=1e-3)
    @test isapprox(result.z, -1.5485089774287586; rtol=1e-3)
    @test result.τ == -0.07552758237689744
    @test isapprox(result.s, -257.5832407750216; rtol=1e-3)
    @test isapprox(result.var_s, 27669.910277723324; rtol=1e-3)
    @test isapprox(result.slope, -0.6382978723404256; rtol=1e-3)
    @test isapprox(result.intercept, 473.29787234042556; rtol=1e-3)
end
