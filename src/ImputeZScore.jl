"""
This module simulates under the generalized linear model (GLM) and
"""
module ImputeZScore

include("sumstats.jl")

export impute_zscore, filter_input!

using DataFrames, SnpArrays

typealias PartitionType Union{Array{Tuple{Int64,Int64},1}, Int64}

"""
Impute association statistics
"""
function impute_zscore(zsc_t::DataFrame, refpanel::SnpData;
  partition::PartitionType = 1000000, buf_size::Int64 = 250000,
  λ::Float64 = 0.1, min_nt::Int64 = 10, r2pred_thres = 0.6)

  # extract legend and genotype matrix
  legend = DataFrame(rsID = refpanel.snpid, pos = refpanel.basepairs,
    A0 = refpanel.allele1, A1 = refpanel.allele2)
  genmat = refpanel.snpmatrix

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
  zsc_out = DataFrame()
  for i = 1:size(partition, 1)

    # get start and end position of the window
    start = max(partition[i][1] - buf_size, min_pos)
    stop = min(partition[i][2] + buf_size, max_pos)

    # load data in the window
    zsc_t_buf = zsc_t[(zsc_t[:pos] .>= start) & (zsc_t[:pos] .< stop), :]
    legend_buf = legend[(legend[:pos].>=start) & (legend[:pos] .< stop), :]
    nsnp, ntyped = size(legend_buf, 1), size(zsc_t_buf, 1)
    nimpute = nsnp - ntyped

    # check if the number of typed snps is too small
    if size(zsc_t_buf, 1) < min_nt continue end

    # find indices of typed snps, snps to impute, and snps to keep
    typed_idx = [legend_buf[:rsID][j] in zsc_t_buf[:rsID] for j=1:nsnp]
    typed_idx = BitArray(typed_idx)
    impute_idx = ~typed_idx

    # compute ld for the locus
    genmat_buf = genmat[:, (legend[:pos] .>= start) & (legend[:pos] .< stop)]
    genmat_buf = convert(Array{Float64,2}, genmat_buf)
    ld_buf = calc_ld(genmat_buf)

    # compute the conditional expectation
    ld_buf_tt = ld_buf[typed_idx, typed_idx] + λ * eye(ntyped)
    ld_buf_tt_inv = pinv(ld_buf_tt)
    ld_buf_it = ld_buf[impute_idx, typed_idx]
    zsc_imp = ld_buf_it * (pinv(ld_buf_tt) * zsc_t_buf[:Z])

    # compute the imputation accuracy
    r2pred = [(ld_buf_it[j,:]'*ld_buf_tt_inv*ld_buf_it[j,:])[1] for j=1:nimpute]
    zsc_imp = zsc_imp ./ sqrt(r2pred)

    # create data frame for imputed z scores
    zsc_imp_buf = DataFrame(rsID   = legend_buf[:rsID][impute_idx],
                            pos    = legend_buf[:pos][impute_idx],
                            A0     = legend_buf[:A0][impute_idx],
                            A1     = legend_buf[:A1][impute_idx],
                            Z      = zsc_imp,
                            r2pred = r2pred)

    # only keep snps in the window
    zsc_imp_buf = zsc_imp_buf[((zsc_imp_buf[:pos] .>= partition[i][1]) &
                               (zsc_imp_buf[:pos] .< partition[i][2])  &
                               (zsc_imp_buf[:r2pred] .>= r2pred_thres)), :]

    # concat results
    zsc_out = [zsc_out; zsc_imp_buf]

  end

  return zsc_out

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
Calculate the LD matrix
"""
function calc_ld(genmat)
  ld = cor(genmat)
  ld[isnan(ld)] = 0.0
  return ld
end

end
