# Input Format

This page describes the format of the input to ImputeZScore.

## GWAS summary association data

For each chromosome, ImputeZScore requires a single text file containing the
following columns in the input GWAS summary association data:

* SNP - rs ID of the SNP (e.g. rs62442).
* BP - Base pair position of the SNP.
* A1 - Effect allele of the SNP. The sign of the Z-score is with respect to this allele.
* A2 - The other allele of the SNP.
* Z - The Z-score of the SNP.

## Genome partition

The genome partition file should be in bed format, one for each chromosome. 

## Reference panel

Reference panels should be in [PLINK format](https://www.cog-genomics.org/plink/2.0/input#bed).

The following is a list of popular publicly available reference panels.

* [1000 Genomes Project](http://www.internationalgenome.org/data/)
* [UK10K](https://www.uk10k.org/data_access.html)

We provide 1000 Genomes reference panel for Europeans [here](https://ucla.box.com/s/4ya4lxvwjsujt6rrotn3keaiw9dkolld).
All SNPs in this reference panel have minor allele frequency greater than 1%.
