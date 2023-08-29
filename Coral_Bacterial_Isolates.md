# CORAL BACTERIAL ISOLATE ANALYSIS
# OVERVIEW
1. Download Download Coral-Associated Bacterial Genomes from NCBI & IMG JGI
2. Run CheckM on Coral Isolate Genomes
3. Run Vibrant on Coral Isolate Genomes
4. Dereplicate Viruses with Virathon
 
# 1. Download Coral-Associated Bacterial Genomes
```bash
# Apparently I have no notes on how this was done
```
# 2. Run CheckM on Coral Isolate Genomes
```bash
# Setup CheckM:
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvzf checkm_data_2015_01_16.tar.gz
checkm data setRoot /path/to/checkmdb/
```
```bash
# Run CheckM
genomes=/path/to/genomes
output=/path/to/output/directory
checkm lineage_wf $genomes $output -t 40 --pplacer_threads 40
```
# 3. Run Vibrant on Coral Isolate Genomes
```bash
# Note: prep files for Vibrant by replacing " " and ":" with "_" and by linearizing the fasta files

# Run Vibrant
input=/path/to/clean/fastas
output=/path/to/output/directory
for f in $input/*.fasta; do
  name=$(basename $f .fasta)
  VIBRANT_run.py -i $f -t 64 -folder $output/${name}_vibrant
done
```
# 4. Dereplicate Viruses with Virathon
```bash
# Note: first concatenate all viruses identified by Vibrant into a single .fna file

# Run Virathon
viruses=/path/to/combined/phages/ALL_CRL_ISOLATE_VIRUSES.fna
python3 Virathon.py --genome_files $viruses --make_pops True --threads 24
```
