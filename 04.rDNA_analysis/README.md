
# Extract ONT UL reads containing rDNA sequences

export PATH=~/miniconda3/envs/ribotinv1.3/bin:$PATH

ribotin-ref --approx-morphsize 45000 -r chm13_rDNAs.fa -t 24 -i haplotype-Mat/Pat_ONTUL.fasta.gz -o res_mat/pat_ont_rDNA

# Map extracted ONT UL reads to the major morphs of CHM13-T2T

minimap2 -ax map-ont -H -k 19 -B 5 -t 48 chm13_major_morphs.fa(https://github.com/maickrau/ribotin_paper_experiments) res_mat/pat_ont_rDNA.fa | samtools view -b -o mat/pat_ont.bam

samtools sort -@ 48 -o mat/pat_ont_sort.bam mat/pat_ont.bam

samtools index -@ 48 mat/pat_ont_sort.bam

# Calculate ONT UL read counts and base counts of each major morphs of CHM13-T2T

bash Read_Base_each_morph_Count.sh

# Calculate the rDNA copy numbler of each major morphs of CHM13-T2T

ddPCR copy numbler * (morph-specific ONT reads base / total ONT base of rDNA)

# Gaps fill in the genome assembly

python rDNA_CN_fill_genome_assembly_v2.py -i1 chr13/14/15/21/22_genome_assembly.fa -i2 chm13_chr13_morph.fa -b rDNA_start-rDNA_end(predicted by barrnap v0.9) -o chr13/14/15/21/22_genome_assembly_T2T.fa -fold rDNA_copy_numbler
