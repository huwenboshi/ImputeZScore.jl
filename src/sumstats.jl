using DataFrames

# define equivalent alleles
equiv = Dict{String, Set{String}}()
equiv["AC"] = Set{String}(["TG", "AC", "TC", "AG"])
equiv["AG"] = Set{String}(["TC", "AG", "TG", "AC"])
equiv["CA"] = Set{String}(["GT", "CA", "GA", "CT"])
equiv["CT"] = Set{String}(["GA", "CT", "GT", "CA"])
equiv["TC"] = Set{String}(["AG", "TC", "AC", "TG"])
equiv["TG"] = Set{String}(["AC", "TG", "AG", "TC"])
equiv["GA"] = Set{String}(["CT", "GA", "CA", "GT"])
equiv["GT"] = Set{String}(["CA", "GT", "CT", "GA"])

# define reversed alleles
reverse = Dict{String, Set{String}}()
reverse["AC"] = Set{String}(["GT", "CA", "CT", "GA"])
reverse["AG"] = Set{String}(["CT", "GA", "GT", "CA"])
reverse["CA"] = Set{String}(["TG", "AC", "AG", "TC"])
reverse["CT"] = Set{String}(["AG", "TC", "TG", "AC"])
reverse["TC"] = Set{String}(["GA", "CT", "CA", "GT"])
reverse["TG"] = Set{String}(["CA", "GT", "GA", "CT"])
reverse["GA"] = Set{String}(["TC", "AG", "AC", "TG"])
reverse["GT"] = Set{String}(["AC", "TG", "TC", "AG"])

# define strand ambiguous alleles
ambiguous = Set{String}(["AT", "CG", "TA", "GC"])

# filter out snp
function filter_input!(zsc_t::DataFrame, legend::DataFrame)

    # create a dictionary of snpid -> ref alt alleles for reference panel
    refpanel_alleles = Dict{String, String}()
    for i = 1:size(legend, 1)
        snpid = legend[:rsID][i]
        alleles = legend[:A0][i] * legend[:A1][i]
        refpanel_alleles[snpid] = alleles
    end

    # count number of entries with the same rs id
    snpid_cnt = Dict{String, Int64}()
    for i = 1:size(zsc_t, 1)
        snpid = zsc_t[:rsID][i]
        if !(snpid in keys(snpid_cnt))
            snpid_cnt[snpid] = 1
        else
            snpid_cnt[snpid] += 1
        end
    end

    # count number of entries with same pos
    pos_cnt = Dict{Int64, Int64}()
    for i = 1:size(zsc_t, 1)
        pos = zsc_t[:pos][i]
        if !(pos in keys(pos_cnt))
            pos_cnt[pos] = 1
        else
            pos_cnt[pos] += 1
        end
    end

    # find snps to remove
    filter_set = Set{Int64}()
    rev_set = Set{Int64}()
    for i = 1:size(zsc_t, 1)

        # extract snp info
        snpid = zsc_t[:rsID][i]
        pos = zsc_t[:pos][i]
        alleles = zsc_t[:A0][i] * zsc_t[:A1][i]

        # check if snp is in reference panel
        if snpid in keys(refpanel_alleles)
            ref_alleles = refpanel_alleles[snpid]

            # snps with matching alleles, do nothing
            if ref_alleles in keys(equiv) && in(alleles, equiv[ref_alleles])
                # do nothing
            # snps with inverse alleles, multiply z-score with -1
            elseif ref_alleles in keys(reverse) &&
                   in(alleles, reverse[ref_alleles])
                push!(rev_set, i)
            # strnad ambiguous snps and snps with no-matching alleles
            else
                push!(filter_set, i)
            end
        # snp not found in reference panel
        else
            push!(filter_set, i)
        end

        # check for snp with duplicate snp id
        if snpid_cnt[snpid] > 1
            push!(filter_set, i)
        end

        # check snp with duplicate position
        if pos_cnt[pos] > 1
            push!(filter_set, i)
        end
    end

    # filter out rows in filter set
    rev_set = sort(collect(IntSet(rev_set)))
    zsc_t[:Z][rev_set] = -1.0*zsc_t[:Z][rev_set]
    filter_set = sort(collect(IntSet(filter_set)))
    deleterows!(zsc_t, filter_set)

    # make the position match
    for i=1:size(zsc_t, 1)
        snpid = zsc_t[:rsID][i]
        if snpid in legend[:rsID]
            zsc_t[:pos][i] = legend[:pos][legend[:rsID] .== snpid][1]
        end
    end
end
