# Impute Z Score

## Getting Started

## Background

In a traditional GWAS, one collects genotype data at a small subset of SNPs
over some individuals, then imputes genotypes across the entire genome, and
finally computes association statistics (e.g. Z-scores) for each genotyped
and imputed SNPs. This procedure can take tremendous amount of time as
genotype imputation is computationally extensive.

Here, we impute Z-scores of ungenotyped SNPs directly from Z-scores of
genotyped SNPs, without first performing genotype imputation, saving hundreds
of hours of CPU time. The idea behind this approach is that Z-scores of
genotyped and ungenotyped SNPs follow a multivariate normal distribution with
LD matrix, which can be estimated from a reference panel, as the covariance
structure -- one can impute the Z-scores of ungenotyped SNPs as the
expectation of Z-scores of ungenotyped SNPs conditional on the Z-scores of
genotyped SNPs.

In detail, let \\(Z = (Z_t, Z_u)\\) be the Z-score vector partitioned into two
components, genotyped (\\(Z_t\\)) and ungenotyped (\\(Z_u\\)). It has been
previously shown that \\(Z\\) has the following distribution,
$$
\left[
\begin{array}{c}
Z_t \\
Z_u
\end{array}
\right]
\sim
MVN
\left(
\left[
\begin{array}{c}
\Lambda_t \\
\Lambda_u
\end{array}
\right]
,
\left[
\begin{array}{cc}
\Sigma_{tt} & \Sigma_{tu}\\
\Sigma_{ut} & \Sigma_{uu}
\end{array}
\right]
\right),
$$
where \\(\Lambda = (\Lambda_t, \Lambda_u)\\) is the non-centrality parameter,
\\(\Sigma_{tt}\\) the LD between genotyped SNPs, \\(\Sigma_{tu}\\) the LD
between genotyped and ungenotyped SNPs, \\(\Sigma_{uu}\\) the LD between
ungenotyped SNPs.

The conditional expectation of \\(Z_u\\) given \\(Z_t\\) is then
$$
Z_u | Z_t
\sim
MVN
\left(
\Lambda_u + \Sigma_{ut} \Sigma^{-1}_{tt} Z_t
,
\Sigma_{uu} - \Sigma_{ut} \Sigma^{-1}_{tt} \Sigma_{tu}
\right).
$$

We impute the Z-scores of ungenotyped SNPs as
\\(\hat{Z}_u = E[Z_u | Z_t] = \Sigma_{ut} \Sigma^{-1}_{tt} Z_t\\) under the
null assumption that \\(\Lambda_u = 0\\). Let
\\(W = \Sigma_{ut} \Sigma^{-1}_{tt}\\). This can be viewed as the weights on
Z-scores of genotyped SNPs in the imputation of Z-scores of ungenotyped SNPs.
Then 
$$
\hat{Z}_u \sim MVN (0, \Sigma_{ut} \Sigma^{-1}_{tt} \Sigma_{tu}),
$$
where each entry $\hat{Z}_{u,i}$ of $\hat{Z}_u$ follows
$$
\hat{Z}_{u,i} \sim N (0, \Sigma_{ut,i*} \Sigma^{-1}_{tt} \Sigma_{tu,*i}).
$$
Here, \\(\Sigma_{ut,i*}\\) denotes the \\(i\\)-th row of \\(\Sigma_{ut}\\) and
\\(\Sigma_{tu,*i}\\) the \\(i\\)-th column of \\(\Sigma_{tu}\\). To obtain an
associations statistics that has mean 0 and variance 1, we standardize
\\(\hat{Z}_{u,i}\\) by
\\(\sqrt{\Sigma_{ut,i*} \Sigma^{-1}_{tt} \Sigma_{tu,*i}}\\). More
specifically, the final imputated association statistics of each SNP is
$$
\hat{Z}_{imp,i} = {\hat{Z}_{u,i} \over \sqrt{\Sigma_{ut,i*} \Sigma^{-1}_{tt} \Sigma_{tu,*i}}}
\sim N(0, 1).
$$

In practice, inverting a large matrix can be time-consuming. Instead, we adopt
a window-based approach, i.e. we impute Z-scores of ungenotyped SNPs one
window at a time.
