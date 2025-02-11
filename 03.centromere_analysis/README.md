
# Predicting the position of the centromeric regions using CenMap

export PATH=~/software/RepeatMasker:$PATH
export PATH=~/miniconda3/envs/R4.1.1/bin:$PATH
export PATH=~/miniconda3/envs/rustybam/bin:$PATH
export PATH=~/software/dna-nn:$PATH
export PATH=~/miniconda3/envs/cenmap_dev/bin:$PATH

snakemake -c 64 -p --workflow-profile none --configfile genome/config/config_CQ.yaml --show-failed-logs --conda-cleanup-pkgs cache

# StainedGlass plots of centromeric regions

export PATH=~/miniconda3/envs/snakemake8.25.5/bin:$PATH
export PATH=~/miniconda3/envs/R/bin:$PATH
export PATH=~/miniconda3/envs/python_and_cli_env/bin:$PATH

snakemake --cores 48 make_figures --configfile=config/config.yaml

# Generates monomeric annotations of HOR of human alpha satellites using HumAS-HMMER_for_AnVIL and StV

nohup bash hmmer-run.sh centromereFasta_split/ AS-HORs-hmmer3.3.2-120124.hmm 60 &
bash stv.sh outputFile_from_HumAS.bed

# Calculate HOR length

export PATH=~/miniconda3/envs/cenmap_dev/bin:$PATH

censtats length --input stv_mat_raw.bed --bp_jump_thr 100000 --arr_len_thr 30000 > stv_mat_HORlength.txt

# Plot HOR length

Rscript HOR_length_plot.R --input GPB_102samples-HORlength.txt --input_A stv_CQv3.0_mat_HORlength.txt --input_B stv_CQv3.0_pat_HORlength.txt

# Plot stv HOR annotations

Rscript plot_cens_stvHOR_update_v3.R -i chm13.stv_raw.bed --input1 CQ_mat.stv_raw.bed --input2 CQ_pat.stv_raw.bed -c chrX -o plot_stvHOR_chm13-cqmat-cqpat_sort_resize.pdf --plot_width 5 --plot_height 12

