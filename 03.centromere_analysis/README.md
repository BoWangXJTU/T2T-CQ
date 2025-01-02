
# Predicting the position of the centromeric regions using CenMap
export PATH=~/software/RepeatMasker:$PATH
export PATH=~/miniconda3/envs/R4.1.1/bin:$PATH
export PATH=~/miniconda3/envs/rustybam/bin:$PATH
export PATH=~/software/dna-nn:$PATH
export PATH=~/miniconda3/envs/cenmap_dev/bin:$PATH

snakemake -c 64 -p --workflow-profile none --configfile genome/config/config_CQ.yaml --show-failed-logs --conda-cleanup-pkgs cache

# Generates monomeric annotations of HOR of human alpha satellites using HumAS-HMMER_for_AnVIL and StV

nohup bash hmmer-run.sh centromereFasta_split/ AS-HORs-hmmer3.3.2-120124.hmm 60 &
bash stv.sh outputFile_from_HumAS.bed

