
###verkko trio###
~/software/meryl/meryl-1.4.1/bin/meryl count compress k=21 memory=200 threads=16 output pat-compress.k21mer.meryl pat_R*.fastq.gz

~/software/meryl/meryl-1.4.1/bin/meryl count compress k=21 memory=200 threads=16 output mat-compress.k21mer.meryl mat_R*.fastq.gz

~/software/meryl/meryl-1.4.1/bin/meryl count compress k=21 memory=200 threads=16 output child-compress.k21mer.meryl child_R*.fastq.gz

bash ~/software/merqury/trio/hapmers.sh mat-compress.k21mer.meryl/ pat-compress.k21mer.meryl/ child-compress.k21mer.meryl/

verkko -d trio --hifi child_HiFi_Reads.gz --nano child_ONTUL_reads.gz --hap-kmers mat-compress.k21mer.hapmer.meryl pat-compress.k21mer.hapmer.meryl trio --threads 144 --local --local-memory 2000 --min-ont-length 100000

###hifiasm trio###
pat_list=`cat paternal.HiFi.fofn`

mat_list=`cat maternal.HiFi.fofn`

yak count -b37 -t48 -o pat.yak $pat_list

yak count -b37 -t48 -o mat.yak $mat_lis

hifiasm -o child.asm --write-ec -t 144 --ul child_ONTUL_100kb_reads.gz -i mat.yak -2 pat.yak child_HiFi_Reads.gz

###hifiasm ONT mode###
hifiasm -o mat -t144 -l0 --ont binned_ONTUL_reads.fastq.gz

###Gap filling###
tgsgapcloser --scaff genome.fasta --reads "binned_ONTUL_reads_100kb.fastq.gz" or "assemblies from hifiasm with ONT mode" --output res --ne --minmap_arg '-x asm20' --thread 60

###Check genome assembly status###
python genome_statistics.py -i genome_final.fasta -o stats.txt

###Juicer Hi-C###
genome=genome_final.fasta

bwa index $genome

python generate_site_positions.py DpnII name $genome

bash ~/software/juicer-1.6/scripts/juicer.sh -t 42 -g name -d name -s "DpnII" -p genome.size -y name_DpnII.txt -z $genome -D ~/software/juicer-1.6

