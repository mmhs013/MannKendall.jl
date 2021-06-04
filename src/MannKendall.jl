module MannKendall

using Statistics
using StatsFuns

export
    mk_original_test

include("utils.jl")
include("MK_tests.jl")

end
