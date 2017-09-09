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
     partition::PartitionType = 1000000, buf_size::Int64 = 250000,
     λ::Float64 = 0.1, min_nt::Int64 = 10)

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

     # iterate through the partition
     for i = 1:size(partition, 1)

          # get start and end position of the window
          start = partition[i][1]
          stop = partition[i][2]

          # add buffer to the window
          start_buf = max(start - buf_size, min_pos)
          stop_buf = min(stop + buf_size, max_pos)

          # load data in the window
          zsc_t_sub = zsc_t[(zsc_t[:pos] .>= start) & (zsc_t[:pos] .< stop), :]
          zsc_t_buf = zsc_t[((zsc_t[:pos] .>= start_buf) &
                             (zsc_t[:pos] .< stop_buf)), :]

          # check if the number of typed snps is too small
          if size(zsc_t_buf, 1) < min_nt
               continue
          end

          # extract reference panel
          legend_sub = legend[(legend[:pos].>=start) & (legend[:pos].<stop),:]
          legend_buf = legend[((legend[:pos].>=start_buf) &
                               (legend[:pos].<stop_buf)),:]
          genmat_buf = convert(Array{Float64,2}, refpanel.snpmatrix[:,
               ((legend[:pos].>=start_buf) & (legend[:pos].<stop_buf))])
          nsnp_all = size(genmat_buf, 1)
          nsnp_typed = size(zsc_t_buf, 1)

          # get untyped snp id
          snpid_untyped = setdiff(Set{String}(legend_buf[:rsID]),
               zsc_t_buf[:rsID])
          snpid_untyped = collect(snpid_untyped)


          # get the array of snp id
          snpid_all = [ Symbol("$snpid") for snpid in legend_buf[:rsID] ]
          snpid_typed = [ Symbol("$snpid") for snpid in zsc_t_buf[:rsID] ]

          # calculate the ld
          ld_buf = calc_ld(genmat_buf)
          ld_buf = convert(DataFrame, ld_buf)
          names!(ld_buf, snpid_all)

          # get the typed ld
          

          break

          #println(ld_buf)

     end
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

"""
Calculate the LD matrix
"""
function calc_ld(genmat)
     ld = cor(genmat)
     ld[isnan(ld)] = 0.0
     return ld
end

end
