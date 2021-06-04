using MannKendall
using Test
# NoTrendData = ones(360)*rand(Int, 1)
TrendData = 1:100;
@testset "MannKendall.jl" begin
    # Write your tests here.
    TrendRes = mk_original_test(TrendData);

    @test TrendRes.trend == "increasing";
    @test TrendRes.h == true;
    @test TrendRes.p == 0.0;
    @test TrendRes.z == 14.73869998208331;
    @test TrendRes.τ == 1.0;
    @test TrendRes.s == 4950;
    @test TrendRes.var_s ==112750.0;

end
