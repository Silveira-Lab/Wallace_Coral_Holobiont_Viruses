# CORAL METAGENOME ANALYSIS
# OVERVIEW
1. Download Coral Metagenomes -  <i>SRA-Toolkit</i>
2. Adapter Trimming, Quality Filtering, Entropy Filtering - <i>BBDuk</i>
3. Classify Reads - <i>Kaiju</i>
4. Assemble Reads - <i>SPAdes</i>
5. Map Reads to Contigs - <i>Bowtie 2</i>
6. Bin bMAGs - <i>CONCOCT, MetaBAT 2</i>
7. Refine bMAGs - <i>CONCOCT, MetaBAT 2</i>
8. Obtain Viral Contigs - <i>Vibrant</i>
9. Bin vMAGs - <i>vRhyme</i>
10. Dereplicate Viruses - <i>BlastN, Virathon</i>
11. Virus Sequence Similarity Tree - [ITOL]
12. Virus Taxonomy - <i>vConTACT2.0</i>
 
# 1. Download Coral Metagenomes
```bash
# Download based on sra accession number [sratoolkit.3.0.2]: 
prefetch --option-file {sra_accession_list}.txt --max-size 10t
```
```bash
# Transform .sra files to .fastq files [sratoolkit.3.0.2]: 
fasterq-dump --split-files {sra_accession}.sra
```

# 2. Adapter Trimming, Quality Filtering, & Entropy Filtering
```bash
# Repair reads: 
repair.sh in1=${name}_R1.fastq in2=${name}_R2.fastq out1=${name}_fixed_R1.fastq out2=${name}_fixed_R2.fastq outs=${name}.singletons.fastq repair
```
```bash
# BBDuk QC - Paired Reads: 
for f in *_fixed_R1.fastq; do
  name=$(basename $f _fixed_R1.fastq)
  bbduk.sh -Xmx512m -da \
  in1=$wd/${name}_fixed_R1.fastq in2=$wd/${name}_fixed_R2.fastq \
  out1=${name}_1out1.fastq out2=${name}_2out1.fastq \
  ktrim=rl k=23 mink=11 hdist=1 qtrim=rl ref=/path/to/adapters.fa qtrim=rl trimq=30 maq=30 entropy=0.90
done

# BBDuk QC - Single Reads: 
for f in *_fixed_R1.fastq; do
  name=$(basename $f _fixed_R1.fastq)
  bbduk.sh -Xmx512m -da in=${name}_fixed_R1.fastq \
  out=${name}_1out1.fastq \
  ktrim=rl k=23 mink=11 hdist=1 qtrim=rl ref=/path/to/adapters.fa qtrim=rl trimq=30 maq=30 entropy=0.90
done
```
```bash
# Visualize QC'd Reads - [FastQC, MultiQC]
for f in *_1out1.fastq; do
  name=$(basename $f _1out1.fastq)
  fastqc {name}.fastq
done

multiqc /path/to/fastqc/directory
```
# 3. Classify Reads
```bash
# Kaiju read taxonomy - Paired reads
db=/path/to/kaiju-db/progenomes_100922/kaiju_db_progenomes.fmi
nodes=/path/to/kaiju-db/progenomes_100922/nodes.dmp
for R1 in $reads/*_1out1.fastq; do
  name=$(basename $R1 _1out1.fastq)
  kaiju -t $nodes -f $db -i ${name}_1out1.fastq -j ${name}_2out1.fastq -o ${name}_kaiju.out
done

# Kaiju read taxonomy - Single reads
db=/path/to/kaiju-db/progenomes_100922/kaiju_db_progenomes.fmi
nodes=/path/to/kaiju-db/progenomes_100922/nodes.dmp
for R1 in $reads/*_1out1.fastq; do
  name=$(basename $R1 _1out1.fastq)
  kaiju -t $nodes -f $db -i ${name}_1out1.fastq -o ${name}_kaiju.out
done
```
```bash
# Kaiju 2 table
input=/path/to/kaiju/output
t=/path/to/kaiju-db/progenomes_100922/nodes.dmp
n=/path/to/kaiju-db/progenomes_100922/names.dmp
for f in $input/*_kaiju.out; do
  name=$(basename $f _kaiju.out)
  kaiju2table -t $t -n $n -u -r genus -l superkingdom,phylum,class,order,family,genus,species -o ${name}_reads_phyla.tsv ${name}_kaiju.out
done
```
# 4. Assemble Reads
```bash
# SPAdes Assembly - Paired Reads: 
for f in *_1out1.fastq; do
  name=$(basename $f _1out1.fastq)
  spades.py --meta --pe1-1 ${name}_1out1.fastq  --pe1-2 ${name}_2out1.fastq --only-assembler -m 100 -o ${name}.spades-out
done
       
# SPAdes Assembly - Single Reads:
for f in *_1out1.fastq; do
  name=$(basename $f _1out1.fastq)
  spades.py -s ${name}_1out1.fastq --only-assembler -o ${name}.spades-out
done
```
# 5. Map Reads to Contigs
```bash
# Create Bowtie2 DB for each sample's contigs
contigs=/path/to/contigs
output=/path/to/output/directory
for f in $contigs/*.contigs.fasta.gz; do
  name=$(basename $f .contigs.fasta.gz)
  bowtie2-build $f $output/${name}_BOWTIE2_CONTIGS_DB
done
```
```bash
# Map quality controlled reads to their own contigs database
reads=/path/to/reads
output=/path/to/output/directory
for f in $reads/*_R1.fastq.gz ; do
  name=$(basename $f _R1.fastq.gz)
  db=/projectnb/viralecology/projects/Wallace/bMAGs/mapping/contigs_dbs/${name}_BOWTIE2_CONTIGS_DB
  bowtie2 --mp 4 -X 1000 -x $db -1 $reads/${name}_R1.fastq.gz  -2 $reads/${name}_R2.fastq.gz  -S $output/${name}.sam
done
```
```bash
# Convert .sam outputs to .sorted.bam 
sams=/path/to/bowtie/sam/files
bams=/path/to/output/bam/directory
for f in $sams/*.sam ; do
  name=$(basename $f .sam)
  samtools sort -o $bams/${name}.align.sort.bam $sams/${name}.sam
  samtools index ${name}.align.sort.bam
done
```
# 6. Bin Bacterial MAGs (bMAGs)
<font size=”6”> CONCOCT </font>
```bash
# Split the reads from your contgs into 10KB chunks and create a bed file
input=/path/to/contigs/fastas
output=/path/to/output/directory
for f in $input/*.contigs.fasta.gz; do
  name=$(basename $f .contigs.fasta.gz)
  cut_up_fasta.py $f -c 10000 -o 0 --merge_last -b $output/${name}_contigs_10K.bed > $output/${name}_contigs_10K.fa
done
```
```bash
# Generate the coverage table
input1=/path/to/bed/files
input2=/path/to/sorted/bam/files
output=/path/to/output/directory
for f in $input1/*_contigs_10K.bed; do
  name=$(basename $f _contigs_10K.bed)
  concoct_coverage_table.py $f $input2/${name}.aln.sort.bam > $output/${name}_coverage_table.tsv
done
```
```bash
# Run CONCOCT
```
```bash
# Merge subcontig clustering into original contig clustering
```
```bash
# Extract bins as individual FASTAs
```
<font size=”6”> MetaBAT 2 </font>
```bash
#ENTER CODE HERE
```
# 7. Refine Bacterial MAGs (bMAGs)
```bash
#ENTER CODE HERE
```
# 8. Obtain Viral Contigs
```bash
# Vibrant Documentation: https://github.com/AnantharamanLab/VIBRANT
```
```bash
# Vibrant viral contig identification
input=/path/to/contigs
output=/path/to/output/directory
for f in $input/*.contigs.fasta; do
name=$(basename $f .contigs.fasta)
VIBRANT_run.py -i $f -t 64 -folder $output/${name}_vibrant
done
```
# 9. Bin Viral MAGs (vMAGs)
```bash
# First, reformat names of forwared and reverse read files
# Forward: _1.fastq or _R1.fastq
# Reverse: _2.fastq or _R2.fastq 

wd=/path/to/viral/contigs
for f in $wd/*_1out1.fastq; do
name=$(basename $f _1out1.fastq)
mv ${name}_1out1.fastq ${name}_R1.fastq
mv ${name}_2out1.fastq ${name}_R2.fastq
done
```
```bash
input=/path/to/viral/contigs/{sample}.fna
output=/path/to/output/directory
reads=/path/to/QC/{sample}_reads
vRhyme -i $input -o $output -u $reads/*.fastq --method longest
conda deactivate
```
# 10. Dereplicate Viruses
<font size=”6”> Dereplicate contigs within bins with BlastN </font> 
```bash
# Make the blast database for your vMAGs
output=/path/to/output/directory
wd=/path/to/vRhyme_best_bins_fasta
for f in $wd/*.fasta; do
	name=$(basename $f .fasta)
	makeblastdb -in ${f} -dbtype nucl -out $out/${name}_db
done
```
```bash
# Run BlastN
wd=/path/to/vRhyme_best_bins_fasta
db=/path/to/vRhyme_best_bins_fasta_database
out=/path/to/output/directory
for f in $wd/*.fasta; do
	name=$(basename $f .fasta)
	blastn -db $db/${name}_db -query $wd/${name}.fasta -out $out/${name}.out -outfmt "6 qseqid sseqid pident qcovs qlen slen length mismatch qstart qend sstart send evalue" -num_threads 16 
done
```
```bash
# Make list of sequence IDs that need to be removed
wd=/path/to/blast/output
for f in $wd/*.out; do
  name=$(basename $f .out)
  awk '{ if ($3>97) print }' $wd/${name}.out | awk '{ if ($1!=$2) print}' | awk '{ if ($5<$6) print}' | awk '{ if ($4>95) print }' | awk '{print $1}' | sort | uniq > $wd/${name}.remove.list
done
```
```bash
# Pull contigs out of the vMAG bins with seqkit
wd=/path/to/blast/output
bins=/path/to/vRhyme_best_bins_fasta
for f in $wd/*.out; do
  name=$(basename $f .out)
  seqkit grep -v -f $wd/${name}.remove.list $bins/${name}.fasta > $bins/no_dups/${name}.no_dups.fasta
done
```
<font size=”6”> Prepare bins for further dereplication </font>
```bash
# Linearize files - convert the fasta files into single line files to N-link
wd=/path/to/dereplicated/bins
for f in $wd/*.no_dups.fasta ; do
  name=$(basename $f .no_dups.fasta)
  awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $wd/${name}.no_dups.fasta >       $wd/single_line_bins_fastas/${name}_single_line.fasta
done
```
```bash
# N-Link contigs within bins
input=/path/to/dereplicated_single_line_bins_fastas
output=/path/to/output/directory
link_bin_sequences.py -i $input -o $output -e fasta -n 1000 -c N
```
```bash
# Make list of contigs that were binned into vMAGs
wd=/path/to/vRhyme/output/directories
for dir in $wd/*; do
	name=$(basename $dir )
	cd $dir
	awk '(NR>1)' $dir/*.membership.tsv | awk '{print $1}' >> ${name}_binned_vContigs.list
	wc -l  ${name}_binned_vContigs.list
done
```
```bash
# Pull binned contigs out to keep unbinned vContigs for database
list=/path/to/list/of/binned/vContigs
fasta=/path/to/concatenated/fasta/file/of/all/viruses
out=/path/to/output/directory
for f in $list/*_binned_vContigs.list; do
	name=$(basename $f _binned_vContigs.list)
	seqkit grep -v -f $list/${name}_binned_vContigs.list $fasta/${name}_ALL_VIRUSES_final.fna > $out/${name}_ALL_VIRUSES_mag-contigs_removed.fna
done
```
<font size=”6”> Dereplicate vMAGs with Virathon </font> 
```bash
# Concatenate N-linked vMAGs
cat *.fasta >> ALL_vMAGS.fasta
```
```bash
# Run Virathon
vMAGs=/path/to/ALL_vMAGS.fasta
Virathon.py --genome_files $vMAGs --make_pops True --threads 24
```
