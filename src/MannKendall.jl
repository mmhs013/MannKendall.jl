module MannKendall

using Statistics
using StatsFuns

include("utils.jl")
include("MK_tests.jl")

# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------
# Sen's slope estimators
export
    sens_slope,
    seasonal_sens_slope,

# Mann-Kendall trend tests
    mk_original_test,
    mk_hamed_rao_modification_test,
    mk_yue_wang_modification_test,
    mk_pre_whitening_modification_test,
    mk_trend_free_pre_whitening_modification_test,
    mk_multivariate_test,
    mk_seasonal_test,
    mk_regional_test,
    mk_correlated_multivariate_test,
    mk_correlated_seasonal_test,
    mk_partial_test

end
