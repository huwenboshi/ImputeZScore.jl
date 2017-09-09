"""
This module simulates under the generalized linear model (GLM) and
"""
module ImputeZScore

include("sumstats.jl")

export impute_zscore

using DataFrames, SnpArrays

typealias PartitionType Union{Array{Tuple{Int64,Int64},1}, Int64}

"""
Impute association statistics
"""
function impute_zscore(zsc_t::DataFrame, refpanel::SnpData;
     partition::PartitionType = 1000000, buf_size::Int64 = 500000,
     maf_th = 0.01, λ::Float64 = 0.1, min_nt::Int64 = 50)

     # create the legend from reference panel
     legend = create_legend(refpanel)

     # create partition if the user specifies window size instead of
     # an array of start and stop tuples
     if typeof(partition) == Int64
          min_pos = minimum(refpanel.basepairs)
          max_pos = maximum(refpanel.basepairs)
          partition = create_partition(min_pos, max_pos, partition)
     end

     # clean up input typed z scores
     filter_input!(zsc_t, legend)

     
end

"""
Create partition based on start and end position and window size
"""
function create_partition(min_pos::Int64, max_pos::Int64, win_sz::Int64)
     nwin = Int64(ceil((max_pos - min_pos + 1) / win_sz))
     windows = Array{Tuple{Int64,Int64},1}(nwin)
     windows[1] = (min_pos, min_pos + win_sz)
     for i = 2:nwin
          start = windows[i-1][2]
          stop = min(windows[i-1][2] + win_sz, max_pos)
          windows[i] = (start, stop)
     end
     return windows
end

"""
Create legend from SnpData
"""
function create_legend(ref_panel::SnpData)
     legend = DataFrame(rsID = ref_panel.snpid, pos = ref_panel.basepairs,
          A0 = ref_panel.allele1, A1 = ref_panel.allele2)
     return legend
end

end
