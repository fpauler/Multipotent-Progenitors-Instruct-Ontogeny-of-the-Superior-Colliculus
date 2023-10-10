# on the cluster environment
# Notes:
mm10_rmsk.gtf: repeatmasker annotation obtained from UCSC genome browser on Nov. 2022
outs: output folder from cellranger analysis

# to get the reference gtf file used by cellranger I did this:
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
tar -zxvf refdata-gex-mm10-2020-A.tar.gz refdata-gex-mm10-2020-A/genes/genes.gtf

# command line to call velocyto:
velocyto run -b ./outs/filtered_feature_bc_matrix/barcodes.tsv.gz -o ./velocyto_out -m mm10_rmsk.gtf ./outs/possorted_genome_bam.bam refdata-gex-mm10-2020-A/genes/genes.gtf

velocyto --version
# velocyto, version 0.17.17

