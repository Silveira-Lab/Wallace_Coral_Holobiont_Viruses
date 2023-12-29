# CORAL METAGENOME ANALYSIS
# OVERVIEW
1. Download Coral Metagenomes -  <i>SRA-Toolkit</i>
2. Adapter Trimming, Quality Filtering, Entropy Filtering - <i>BBDuk</i>
3. Classify Reads - <i>Kaiju</i>
4. Assemble Reads - <i>SPAdes</i>
5. Map Reads to Contigs - <i>Bowtie 2</i>
6. Bin bMAGs - <i>CONCOCT, MaxBin 2, MetaBAT 2</i>
7. Refine & Dereplicate bMAGs - <i>MetaWRAP, anvi'o</i>
8. Obtain Viral Contigs - <i>Vibrant</i>
9. Bin vMAGs - <i>vRhyme</i>
10. Dereplicate Viruses - <i>BlastN, Virathon</i>
11. Assess vMAG Quality - <i>CheckV</i>
12. Virus Sequence Similarity Tree - <i>GL-UVAB</i>
13. Virus Taxonomy - <i>PTT, Kaiju, Kraken</i>
14. Bacterial and Viral Abundances - <i>Smalt</i>
15. Virulence Factor Analysis - <i>BlastP</i>
16. Linking Viruses and their Hosts - <i>MinCED, BlastN, Minimap2, CheckV</i>
 
# 1. Download Coral Metagenomes
```bash
# Download based on sra accession number: 
prefetch --option-file {sra_accession_list}.txt --max-size 10t
```
```bash
# Transform .sra files to .fastq files: 
fasterq-dump --split-files {sra_accession}.sra
```

# 2. Adapter Trimming, Quality Filtering, & Entropy Filtering
```bash
BBDuk Documentation: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/
```
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
Kaiju Documentation: https://github.com/bioinformatics-centre/kaiju
```
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
SPAdes Documentation: https://github.com/ablab/spades
```
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
Bowtie2 Documentation: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
```
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
  samtools sort -o $bams/${name}.align.sort.bam $sams/${name}.sam;
  rm $sams/${name}.sam
  samtools index $bams/${name}.align.sort.bam
done
```
# 6. Bin Bacterial MAGs (bMAGs)
<font size=”6”> CONCOCT </font>
```bash
# CONCOCT Documentation: https://concoct.readthedocs.io/en/latest/index.html#
# Note: don't try to run any of this with zipped files
```
```bash
# Split the reads contgs into 10KB chunks and create a bed file
input=/path/to/contigs/fastas
output=/path/to/output/directory
for f in $input/*.contigs.fasta; do
  name=$(basename $f .contigs.fasta)
  cut_up_fasta.py $f -c 10000 -o 0 --merge_last -b $output/${name}_contigs_10K.bed > $output/${name}_contigs_10K.fa
done
```
```bash
# Generate coverage table
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
comp_file=/path/to/composition/files
cov_file=/path/to/coverage/files
output=/path/to/output/directory
for f in $comp_file/*_contigs_10K.fa; do
  name=$(basename $f _contigs_10K.fa)
  concoct --composition_file $comp_file/${name}_contigs_10K.fa --coverage_file $cov_file/${name}_coverage_table.tsv -b $output/${name}_concoct-bins/ --threads 64
done
```
```bash
# Merge subcontig clustering into original contig clustering
bins=/path/to/bins/directories
for dir in $bins/*_concoct-bins; do
	name=$(basename $dir _concoct-bins)
	cd $dir
	mv clustering_gt1000.csv ${name}_clustering_gt1000.csv
	merge_cutup_clustering.py ${name}_clustering_gt1000.csv > ${name}_clustering_merged.csv
done
```
```bash
# Extract fasta bins
assemblies=/path/to/contigs/fastas
bins=/path/to/concoct/output/bins
out=/path/to/output/directory
for f in $assemblies/*.contigs.fasta; do
	name=$(basename $f .contigs.fasta)
	mkdir $out/${name}_bins
	extract_fasta_bins.py $assemblies/${name}.contigs.fasta $bins/${name}_concoct-bins/${name}_clustering_merged.csv --output_path $out/${name}_bins
done
```
<font size=”6”> MetaBAT 2 </font>
```bash
# Generate Depth Profiles
bams=/path/to/bowtie2/mapping/bams/directory
out=/path/to/output/directory
for f in $bams/*.align.sort.bam; do
	name=$(basename $f .align.sort.bam)
	jgi_summarize_bam_contig_depths --outputDepth $out/${name}_depth.txt $bams/${name}.align.sort.bam
done
```
```bash
# Run MetaBAT 2
contigs=/path/to/contigs/fastas
depths=/path/to/depth_profiles
out=/path/to/output/directory
for f in $contigs/*.fasta; do
	name=$(basename $f .fasta)
	metabat2 -m 1500 -t 64 -i $f -a $depths/${name}_depth.txt -o $out/${name}_metabat2-bins
done
```
<font size=”6”> MaxBin2 </font>
```bash
# Make abundance files from bam mapping outputs (must be tab-delimited format) - (contig header)\t(abundance)
wd=/path/to/mapping/bam/files
for i in $wd/*.align.sort.bam; do
    samtools idxstats $i > $i.idxstats
    cut -f1,3 $i.idxstats > $i.counts.tsv
done
```
```bash
# Run MaxBin2
contigs=/path/to/contigs
abundance=/path/to/abundance/files
out=/path/to/output/directory
for f in $contigs/*.contigs.fasta.gz; do
	name=$(basename $f .contigs.fasta.gz)
	run_MaxBin.pl -thread 64 -contig $f -out $out/${name} -abund $abundance/${name}.align.sort.bam.counts.tsv
done
```
# 7. Refine and Dereplicate Bacterial MAGs (bMAGs)
```bash
# Run metawrap bin refinement https://github.com/bxlab/metaWRAP
metabat_bins=/path/to/metabat/bins
maxbin_bins=/path/to/maxbin/bins
concoct_bins=/path/to/concoct/bins
out=/path/to/output/directory
for f in $concoct_bins/*_bMAGs; do
	name=$(basename $f _bMAGs)
	metawrap bin_refinement -o $out/${name}_bins -t 64 -A $concoct_bins/${name}_bMAGs -B $maxbin_bins/${name}_bMAGs -C $metabat_bins/${name}_bMAGs -c 20 -x 10
done
```
```bash
# Run anvi'o dereplication using a 95% similarity threshold https://anvio.org/help/7/programs/anvi-dereplicate-genomes/
# Note bMAGs_list.tsv is a tab-separated file with two columns. The first column is the name of the bin and the second column is the path to the fasta file.

bins_list=/path/to/bMAGs_list.tsv
out=/path/to/output/directory
anvi-dereplicate-genomes -f $bins_list -o $out --program fastANI --similarity-threshold 0.95
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
```bash
# Blast Documentation: https://www.ncbi.nlm.nih.gov/books/NBK279690/
# Virathon Documentation: https://github.com/felipehcoutinho/virathon
```
<font size=”6”> Dereplicate contigs within vMAG bins with BlastN </font> 
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

# Pull contigs out of the vMAG bins with seqkit
```
<font size=”6”> Prepare vMAG bins for further dereplication </font>
```bash
# Linearize files - convert the fasta files into single line files to N-link
wd=/path/to/dereplicated/bins
for f in $wd/*.no_dups.fasta ; do
  name=$(basename $f .no_dups.fasta)
  awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $wd/${name}.no_dups.fasta >       $wd/single_line_bins_fastas/${name}_single_line.fasta
done
```
```bash
# N-Link contigs within bins using vRhyme's link_bin_sequences.py
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
<font size=”6”> Dereplicate Remaining vContigs </font> 
```bash
fna=/path/to/remaining/vCongigs/ALL_VIRUSES_mag-contigs_removed.fna  
python3 Virathon.py --genome_files $fna --make_pops True --threads 24
```
# 11. Assess vMAG Quality
```bash
# CheckV Documentation: https://bitbucket.org/berkeleylab/checkv/src/master/
```
```bash
# Assess quality of binned vMAGs
fasta=/path/to/concatenated/vMAGs/DEREP_vMAGs.fasta
output=/path/to/output/directory
checkv end_to_end $fasta $output -t 16 -d checkv-db-v1.4
```
# 12. Phylogenomic Tree of Coral Viruses
```bash
# GL-UVAB Documentation: https://sourceforge.net/projects/gluvab/files/
```
<font size=”6”> Select High & Medium Quality Viruses for Tree </font> 
```bash
# Make list of high & medium quality vMAGs using CheckV output & pull that list from original fasta
# Make list of complete circular, high quality, and medium quality vContigs using Vibrant output & pull that list from original fasta
```
<font size=”6”> Dereplicate Viral RefSeq with CheckV's rapid genome clustering based on pairwise ANI </font> 
```bash
# Create Blast Database
out=/path/to/output/directory
wd=/path/to/ViralRefSeq_1.1_2023-07-13
makeblastdb -in $wd/viral.1.1.genomic.fna -dbtype nucl -out viral.1.1.genomic_db
```
```bash
# Use megablast from blast+ package to perform all-vs-all blastn of sequences
out=/path/to/output/directory
wd=/path/to/ViralRefSeq_1.1_2023-07-13
db=/path/to/ViralRefSeq_1.1_2023-07-13/blastdb
blastn -query $wd/viral.1.1.genomic.fna -db $db/viral.1.1.genomic_db -outfmt '6 std qlen slen' -max_target_seqs 10000 -out $out/viral.1.1.genomic_blast.tsv -num_threads 32
```
```bash
# Calculate pairwise ANI by combining local alignments between sequence pairs:
blastout=/path/to/ViralRefSeq_1.1_2023-07-13/MIUVIG_Viral_RefSeq_2023-09-27
anicalc.py -i $blastout/viral.1.1.genomic_blast.tsv -o $blastout/viral.1.1.genomic_ani.tsv
```
```bash
# Perform UCLUST-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF):
wd=/path/to/ViralRefSeq_1.1_2023-07-13
blastout=/path/to/ViralRefSeq_1.1_2023-07-13/MIUVIG_Viral_RefSeq_2023-09-27
aniclust.py --fna $wd/viral.1.1.genomic.fna --ani $blastout/viral.1.1.genomic_ani.tsv --out $blastout/viral.1.1.genomic_clusters.tsv --min_ani 95 --min_tcov 85 --min_qcov 0
```
```bash
# Grab Colum 1 (representative seqs) and make a list:
awk '{print $1}' viral.1.1.genomic_clusters.tsv > MIUVIG_Viral_RefSeq_2023-09-27.tsv
cat MIUVIG_Viral_RefSeq_2023-09-27.tsv | wc -l

# Check for Unique ID's
sort MIUVIG_Viral_RefSeq_2023-09-27.tsv | uniq > MIUVIG_Viral_RefSeq_2023-09-27_unique.tsv
cat MIUVIG_Viral_RefSeq_2023-09-27_unique.tsv | wc -l

# Create New Dereplicated Database
wd=/path/to/ViralRefSeq_1.1_2023-07-13/MIUVIG_Viral_RefSeq_2023-09-27
fasta=/path/to/ViralRefSeq_1.1_2023-07-13/viral.1.1.genomic.fna
seqkit grep -f $wd/MIUVIG_Viral_RefSeq_2023-09-27_unique.tsv $fasta >> $wd/MIUVIG_Viral_RefSeq_2023-09-27_unique.fna
```

<font size=”6”> GL-UVAB Part 1 </font> 
```bash
# GLUVAB Pt 1
gf1=/path/to/dereplicated/ViralRefSeq/database
gf2=/path/to/HQ_MQ_viruses/CHOP_DEREP_HQMQ_viruses.fasta
perl GLUVAB_polyN_v0.6.pl --threads 16 --genomes_file_1 $gf1 --genomes_file_2 $gf2
```
```bash
# Select RefSeq ID's that matched to viruses in my DB
awk '{print $2}' GLUVAB_Scaffold_Recip_Scores.tsv > ALL_RefSeq_forgluvab_unsorted_2.txt
# Make list of unique names
sort ALL_RefSeq_forgluvab_unsorted_2.txt | uniq > ALL_RefSeq_forgluvab_unique_2.txt
# Now grab those IDs from the RefSeq database fasta file
db=/path/to/RefSeq/database/MIUVIG_Viral_RefSeq_2023-09-27_unique.fna
seqkit grep -f ALL_ICTV_forgluvab_unique_2.txt $db > ALL_RefSeq_forgluvab_unique_2.fna

# Combine RefSeq Viruses with My Viral Database
myphages=/path/to/CHOP_DEREP_HQMQ_viruses.fasta
RefSeqphages=/path/to/ALL_RefSeq_forgluvab_unique_2.fna
cat $myphages $ictvphages >> HQMQ_CHoP_x_RefSeq.fasta
```
<font size=”6”> GL-UVAB Part 2 </font>
```bash
# GLUVAB Pt 2
gf1=/path/to/HQMQ_CHoP_x_RefSeq.fasta
perl GLUVAB_polyN_v0.6.pl --threads 16 --genomes_file_1 $gf1
```
```bash
# Visualize tree with ITOL
```

# 13. Virus Taxonomy
<font size=”6”> PTT (Phage Taxonomy Tool) </font>
```bash
# PTT Documentation: https://github.com/AnantharamanLab/Kieft_and_Zhou_et_al._2020
# For this, I created a new PTT_virus_taxonomy.tsv file with updated lineages
```
```bash
# Run PTT
FAA=/path/to/coral/viruses/proteindb/my_viruses.faa
python3 PTT.py -i $FAA -f prot
conda deactivate
```
<font size=”6”> Kaiju </font>
```bash
# Kaiju Documentation: https://github.com/bioinformatics-centre/kaiju
# Using Viral RefSeq for the DB here for updated taxonomy
```
```bash
# Run Kaiju
nodes=/path/to/kaiju-db/viral_refseq_091523/nodes.dmp
fmi=/path/to/kaiju-db/viral_refseq_091523/viruses/kaiju_db_viruses.fmi
input=/path/to/my_db/HQMQ_CHoP_x_RefSeq.fasta
out=/path/to/output/directory
kaiju -t $nodes -f $fmi -i $input -v -z 24 -o $out/kaiju_HQMQ.out
```
```bash
# Run Kaiju-addTaxonNames
nodes=/path/to/kaiju-db/viral_refseq_091523/nodes.dmp
names=/path/to/kaiju-db/viral_refseq_091523/names.dmp
input=/path/to/kaiju/output/kaiju_HQMQ.out
out=/path/to/output/directory
kaiju-addTaxonNames -t $nodes -n $names -i $input -o $out/kaiju_taxonnames_HQMQ.tsv
```
<font size=”6”> Kraken </font>
```bash
# Kraken Documentation: https://github.com/DerrickWood/kraken
```
```bash
# Download Taxonomy
DBNAME=/path/to/ViralRefSeq_1.1_2023-07-13/KRAKEN2_DB/ViralRefSeq_7.13.23
kraken2-build --download-taxonomy --db $DBNAME

# Download Library
DBNAME=/path/to/ViralRefSeq_1.1_2023-07-13/KRAKEN2_DB/ViralRefSeq_7.13.23
kraken2-build  --download-library viral --db $DBNAME

# Build DB
DBNAME=/path/to/ViralRefSeq_1.1_2023-07-13/KRAKEN2_DB/ViralRefSeq_7.13.23
kraken2-build --build --db $DBNAME
```
```bash
# Classify Viruses
db=/path/to/ViralRefSeq_1.1_2023-07-13/KRAKEN2_DB/ViralRefSeq_7.13.23
seqs=/path/to/HQMQ_CHoP_x_RefSeq.fasta
out=/path/to/output/directory
kraken2 --db $db $seqs \
--classified-out $out/classified_viruses.out.fq  \
--unclassified-out $out/unclassified_viruses.out.fq \
--output $out/CHoP_viruses.kraken2.txt --report $out/CHoP_viruses.kraken2_report 
```
```bash
# Get Viral Lineages and Use Tax ID to Pair Viruses with their Linneages assigned by Kraken
# https://github.com/zyxue/ncbitax2lin
# Use Tax Dump to Get Lineages
taxdump=/path/to/taxdump/nodes.dmp
names=/path/to/taxdump/names.dmp
out=/path/to/output/directory/ncbi_lineages.csv.gz
ncbitax2lin --nodes-file $taxdump --names-file $names --output $out

# Next, filter ncbi_lineages.csv.gz output in R and use tax_id from Kraken to match lineages with taxa
```

# 14. Bacterial and Viral Abundances
```bash
Smalt Documentation: https://github.com/rcallahan/smalt
```
<font size=”6”> Fraction of quality controlled reads mapped to viruses and bacteria identified in coral metagenomes </font> 
```bash
# Make Smalt index
out=/path/to/output/directory
db=/path/to/viral/or/bacterial/database/{name}.fasta
smalt index $out/SMALT_INDEX $db

# Run Smalt
reads=/path/to/reads
sm=/path/to/CHoP_SMALT_INDEX
out=/path/to/output/directory
for R1 in $reads/*_R1.fastq; do
	name=$(basename $R1 _R1.fastq)
	R2=${name}_R2.fastq
	smalt map -n 40 -y 0.80 -o $out/${name}.align.sam $sm $R1 $reads/$R2
	samtools sort -o $out/${name}.align.sort.bam $out/${name}.align.sam
	rm $out/${name}.align.bam
	samtools index $out/${name}.align.sort.bam
done
```

# 15. Virulence Factor Analysis
```bash
Blast Documentation: https://rnnh.github.io/bioinfo-notebook/docs/blast.html#the-command-line-version-of-blast
```
```bash
# Download VFDB
wget https://figshare.com/ndownloader/files/15347432
mv 15347432 FileS1_VF_DB.fasta

# Pull out seq names
grep '>' FileS1_VF_DB.fasta > metadata.tsv

# Tidy up the sequence names
sed -e 's/ /_/g' FileS1_VF_DB.fasta | sed -e 's/(//g' | sed -e 's/)//g' | sed -e 's/\[//g' | sed -e 's/\]//g' | sed -e 's/__/_/g' | sed  's/gb|/gb/g' > VFDB.faa
```
```bash
# Make Blast DB out of VFDB
out=/path/to/output/directory
wd=/path/to/VFDB/directory
for f in $wd/VFDB.faa; do
	name=$(basename $f .faa)
	makeblastdb -in ${f} -dbtype prot -out $out/${name}_db
done
```
```bash
# Run BlastP
wd=/path/to/my/DB
db=/path/to/VFDB/blastdb
out=/path/to/output/directory
blastp -query $wd/coral_viruses_protein_DB.faa -db $db/VFDB_db -out $out/coral_VFs_blast.out.tsv -outfmt 6 -evalue 10e-5
```

# 15. Phage Host Linkages
<font size=”6”> CRISPR spacer linkages </font>
```bash
# Run minced on bmags
minced=/path/to/minced
bacteria=/path/to/all/bacterial/genomes/and/MAGs/bacteria.fasta
$minced -spacers $bmags ALL_BACTERIA_CRISPR.txt ALL_BACTERIA_CRISPR.gff
```
```bash
# Make blast DB of Spacers
makeblastdb=/path/to/makeblastdb
wd=/path/to/CRISPR/spacers
out=/path/to/output/directory
$makeblastdb -in $wd/ALL_BACTERIA_CRISPR_spacers.fa -out $out/ALL_BACTERIA_CRISPR_spacers -dbtype 'nucl' -hash_index
```
```bash
# BlastN CRISPR spacers against viral database
blastn=/path/to/blastn
viruses=/path/to/final/viral/database
blastdb=/path/to/CRISPR/spacers/db
out=/path/to/output/directory
$blastn -query $viruses/GCVDB.fasta -db $blastdb/ALL_BACTERIA_CRISPR_spacers -max_target_seqs 1 -task "blastn-short" -outfmt 6 > $out/CRISPR_spacer_links.BLASTn.tsv
```

<font size=”6”> Provirus linkages </font>
```bash
# Map GCVDB viruses to Bacteria using Minimap2
minimap=/path/to/minimap2
bmags=/path/to/bMAGSs.fa
viruses=/path/to/final/viral/database/GCVDB.fasta
out=/path/to/output/directory
$minimap -c $bmags $viruses > $out/minimap_viruses_to_bmag_contigs.paf

# Pull out bacterial contigs that mapped to viruses --> provirus_contigs.fa
```
```bash
# Search for Proviruses with CheckV
checkv=/path/to/checkv
fasta=/path/to/provirus_contigs.fa
output=/path/to/output/directory
db=/path/to/checkv-db-v1.4
$checkv end_to_end $fasta $output -t 16 -d $db

# Check 'contamination.tsv' for proviruses
```



