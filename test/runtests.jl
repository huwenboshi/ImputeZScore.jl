module ImputeZScoreTest

#include("../src/ImputeZScore.jl")

using DataFrames, ImputeZScore, SnpArrays

if VERSION >= v"0.5.0-dev+7720"
    using Base.Test
else
    using BaseTestNext
end

# load test data set

refpanel = SnpData(Pkg.dir("ImputeZScore") * "/docs/1000G.EUR.22");
zsc_t = readtable(Pkg.dir("ImputeZScore") * "/docs/hdl_chr22_typed.zsc",
    separator=' ');
impute_zscore(zsc_t, refpanel)

end
