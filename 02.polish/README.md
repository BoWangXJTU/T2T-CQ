
# classify and bin HiFi reads
yak triobin pat.yak mat.yak child_HiFi_Reads.gz > triobin.out

awk '$2=="p"||$2=="a"||$2==0' triobin.out |cut -f 1 > pat.read.list
fxTools getseq -r pat.read.list child_HiFi_Reads.gz > pat.hifi.fasta

awk '$2=="m"||$2=="a"||$2==0' triobin.out |cut -f 1 > mat.read.list
fxTools getseq -r mat.read.list child_HiFi_Reads.gz > mat.hifi.fasta

# Construct k-mer db
export PATH=~/software/meryl/meryl-1.4.1/bin:$PATH
meryl count k=21 memory=200 threads=48 output child.k21mer.meryl child_R*.fastq.gz

# Collect histogram
meryl histogram child.k21mer.meryl > reads.hist
meryl greater-than 1 child.k21mer.meryl output reads.gt1.meryl

export PATH=~/software/racon/build/bin:$PATH
export PATH=~/miniconda3/envs/winnowmap/bin:$PATH
export PATH=~/miniconda3/envs/merfin/bin:$PATH
export PATH=~/miniconda3/envs/falconc/bin:$PATH
export PATH=~/software/meryl/meryl-1.4.1/bin:$PATH
export PATH=~/miniconda3/envs/yak/bin:$PATH

yak count -t 5 -k 21 -b 37 -o sr.k21.yak child_R*.fastq.gz
yak count -t 5 -k 31 -b 37 -o sr.k31.yak child_R*.fastq.gz

# Polishing with NextPolish2
export PATH=~/miniconda3/envs/nextpolish2/bin:$PATH

bash automated-polishing.sh 36 1 ../draft/cq_t2t_genome_draft.fa pat.hifi.fasta reads.gt1.meryl racon.meryl

samtools index -@ 24 racon.meryl.iter_1.winnowmap.sorted.bam
nextPolish2 -t 48 -r racon.meryl.iter_1.winnowmap.sorted.bam ../draft/cq_t2t_genome_draft.fa sr.k21.yak sr.k31.yak -o ../draft/cq_t2t_genome_draft.np2.fa
