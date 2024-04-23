################
##1. GENOME ASSEMBLY
################
mkdir raw_data
cat /mnt/longship/users/sfo503/Harper_Vault/rice_genomes/X204SC22122095-Z01-F001_01/01.RawData/SRA2_1_C/SRA2_1_C_EKDN230020254-1A_HFCCGDSX7_L3_1.fq.gz X204SC22122095-Z01-F001_01/01.RawData/SRA2_1_C/SRA2_1_C_EKDN230020254-1A_HFCCGDSX7_L4_1.fq.gz >raw_data/SRA2_1_C_1.fq.gz
cat /mnt/longship/users/sfo503/Harper_Vault/rice_genomes/X204SC22122095-Z01-F001_01/01.RawData/SRA2_1_C/SRA2_1_C_EKDN230020254-1A_HFCCGDSX7_L3_2.fq.gz X204SC22122095-Z01-F001_01/01.RawData/SRA2_1_C/SRA2_1_C_EKDN230020254-1A_HFCCGDSX7_L4_2.fq.gz > raw_data/SRA2_1_C_2.fq.gz
cat /mnt/longship/users/sfo503/Harper_Vault/rice_genomes/X204SC22122095-Z01-F001_02/01.RawData/PDP211_B/PDP211_B_EKDN230020255-1A_HFCCGDSX7_L3_1.fq.gz X204SC22122095-Z01-F001_02/01.RawData/PDP211_B/PDP211_B_EKDN230020255-1A_HFCCGDSX7_L4_1.fq.gz > raw_data/PDP211_B_1.fq.gz
cat /mnt/longship/users/sfo503/Harper_Vault/rice_genomes/X204SC22122095-Z01-F001_02/01.RawData/PDP211_B/PDP211_B_EKDN230020255-1A_HFCCGDSX7_L3_2.fq.gz X204SC22122095-Z01-F001_02/01.RawData/PDP211_B/PDP211_B_EKDN230020255-1A_HFCCGDSX7_L4_2.fq.gz > raw_data/PDP211_B_2.fq.gz
cat /mnt/longship/users/sfo503/Harper_Vault/rice_genomes/X204SC22122095-Z01-F001_03/01.RawData/GL37_A/GL37_A_EKDN230020252-1A_HFCCGDSX7_L3_1.fq.gz X204SC22122095-Z01-F001_03/01.RawData/GL37_A/GL37_A_EKDN230020252-1A_HFCCGDSX7_L4_1.fq.gz > raw_data/GL37_A_1.fq.gz
cat /mnt/longship/users/sfo503/Harper_Vault/rice_genomes/X204SC22122095-Z01-F001_03/01.RawData/GL37_A/GL37_A_EKDN230020252-1A_HFCCGDSX7_L3_2.fq.gz X204SC22122095-Z01-F001_03/01.RawData/GL37_A/GL37_A_EKDN230020252-1A_HFCCGDSX7_L4_2.fq.gz > raw_data/GL37_A_2.fq.gz
cat /mnt/longship/users/sfo503/Harper_Vault/rice_genomes/X204SC22122095-Z01-F001_04/01.RawData/J23_A/J23_A_EKDN230020253-1A_HFCCGDSX7_L3_1.fq.gz X204SC22122095-Z01-F001_04/01.RawData/J23_A/J23_A_EKDN230020253-1A_HFCCGDSX7_L3_1.fq.gz > raw_data/J23_A_1.fq.gz
cat /mnt/longship/users/sfo503/Harper_Vault/rice_genomes/X204SC22122095-Z01-F001_04/01.RawData/J23_A/J23_A_EKDN230020253-1A_HFCCGDSX7_L3_2.fq.gz X204SC22122095-Z01-F001_04/01.RawData/J23_A/J23_A_EKDN230020253-1A_HFCCGDSX7_L4_2.fq.gz > raw_data/J23_A_2.fq.gz

#######for each genome do
file=names_raw_data #### genome GL37_A
while read -r line;do
./bbmap/repair.sh -Xmx64G in=/mnt/lustre/groups/bio-ash-2019/rice_genomes/raw_data/{line}_1.fq.gz in2=/mnt/lustre/groups/bio-ash-2019/rice_genomes/raw_data/{line}_2.fq.gz out={line}_sorted.1.fq out2={line}_sorted.2.fq outs={line}_singletons.fq
./bbmap/bbnorm.sh in={line}_sorted.1.fq in2={line}_sorted.2.fq out={line}_normalized1.fq out2={line}_normalized2.fq target=100 
spades.py -k 21,33,55,77 -t 40 --careful -1 {line}_normalized1.fq -2 {line}_normalized2.fq -o {line}_spades_output
done < "file"

mkdir assembly
file=names_raw_data7 #### genome J23_A
while read -r line;do
./bbmap/bbnorm.sh in=/mnt/lustre/groups/bio-ash-2019/rice_genomes/raw_data/{line}_1.fq.gz in2=/mnt/lustre/groups/bio-ash-2019/rice_genomes/raw_data/{line}_2.fq.gz out={line}_sorted.1.fq out2={line}_sorted.2.fq  out=normalized1.fq out=normalized2.fq target=200 
/bbmap/bbnorm.sh in={line}_sorted.1.fq in2={line}_sorted.2.fq out={line}_normalized1.fq out2={line}_normalized2.fq target=100 
spades.py -k 21,33,55,77 -t 40 --careful -1 {line}_normalized1.fq -2 {line}_normalized2.fq -o {line}_spades_output
done < "file"

file=names_raw_data8  #### genome PDP211_B
while read -r line;do
./bbmap/repair.sh -Xmx64G in=/mnt/lustre/groups/bio-ash-2019/rice_genomes/raw_data/{line}_1.fq.gz in2=/mnt/lustre/groups/bio-ash-2019/rice_genomes/raw_data/{line}_2.fq.gz out={line}_sorted.1.fq out2={line}_sorted.2.fq outs={line}_singletons.fq
/bbmap/bbnorm.sh in={line}_sorted.1.fq in2={line}_sorted.2.fq out={line}_normalized1.fq out2={line}_normalized2.fq target=100 
spades.py -k 21,33,55,77 -t 40 --careful -1 {line}_normalized1.fq -2 {line}_normalized2.fq -o {line}_spades_output
done < "file"

file=names_raw_data9 #### genome SRA2_1_C
while read -r line;do
./bbmap/repair.sh -Xmx64G in=/mnt/lustre/groups/bio-ash-2019/rice_genomes/raw_data/{line}_1.fq.gz in2=/mnt/lustre/groups/bio-ash-2019/rice_genomes/raw_data/{line}_2.fq.gz out={line}_sorted.1.fq out2={line}_sorted.2.fq outs={line}_singletons.fq
/bbmap/bbnorm.sh in={line}_sorted.1.fq in2={line}_sorted.2.fq out={line}_normalized1.fq out2={line}_normalized2.fq target=100 
spades.py -k 21,33,55,77 -t 40 --careful -1 {line}_normalized1.fq -2 {line}_normalized2.fq -o {line}_spades_output
done < "file"


###########DO PILON for polishing with the same Illumina data

while read -r line;do
mkdir sibeliaz_{line}
mkdir pilon
bwa index {line}_spades_output/scaffolds.fasta
bwa mem {line}_spades_output/scaffolds.fasta {line}_normalized1.fq {line}_normalized2.fq  -t SLURM_CPUS_PER_TASK > pilon/{line}.sam
samtools view -S -bo -h pilon/{line}.sam -o pilon/{line}.bam
samtools sort pilon/{line}.bam -o pilon/{line}_sorted.bam 
samtools index pilon/{line}_sorted.bam
java -Xmx100G -jar /users/sfo503/scratch/genome/pilon/pilon-1.24.jar --genome  {line}_spades_output/scaffolds.fasta --frags pilon/{line}_sorted.bam --diploid --fix "all" --output pilon/{line}_pilon --changes --threads 30 
seqkit stat -a pilon/{line}_pilon.fasta > pilon/{line}_pilon_sekit.txt
ln pilon/{line}_pilon.fasta ./
done < "file"

####USESIBELIAZ and RACON to MAP against IRGSP_1.0 genome of references (Japonica)
#####modify names in files 
awk '/^>/{print ">genome.seq" ++i; next}{print}' < IRGSP_1.0.genome.fasta > IRGSP_1.0.genome_names.fasta
sed 's/genome/IRGSP/g' IRGSP_1.0.genome_names.fasta > IRGSP_1.0.genome_names2.fasta
sed 's/genome//g' IRGSP_1.0.genome_names.fasta > IRGSP_1.0.genome_names_ragout.fasta

awk '/^>/{print ">genome.seq" ++i; next}{print}' < GL37_A_pilon.fasta > GL37_A_pilon_names.fasta
sed 's/genome/GL37_A/g' GL37_A_pilon_names.fasta > GL37_A_pilon_names2.fasta
sed 's/genome.//g' GL37_A_pilon_names.fasta > GL37_A_pilon_names_ragout.fasta

awk '/^>/{print ">genome.seq" ++i; next}{print}' < PDP211_B_pilon.fasta > PDP211_B_pilon_names.fasta
sed 's/genome/PDP211_B/g' PDP211_B_pilon_names.fasta > PDP211_B_pilon_names2.fasta
sed 's/genome.//g' PDP211_B_pilon_names.fasta > PDP211_B_pilon_names_ragout.fasta

ln pilon/SRA2_1_C_pilon.fasta ./
awk '/^>/{print ">genome.seq" ++i; next}{print}' < SRA2_1_C_pilon.fasta > SRA2_1_C_pilon_names.fasta
sed 's/genome/SRA2_1_C/g' SRA2_1_C_pilon_names.fasta > SRA2_1_C_pilon_names2.fasta
sed 's/genome.//g' SRA2_1_C_pilon_names.fasta > SRA2_1_C_pilon_names_ragout.fasta

awk '/^>/{print ">genome.seq" ++i; next}{print}' < J23_A_pilon.fasta > J23_A_pilon_names.fasta
sed 's/genome/J23_A/g' J23_A_pilon_names.fasta > J23_A_pilon_names2.fasta
sed 's/genome.//g' J23_A_pilon_names.fasta > J23_A_pilon_names_ragout.fasta

sibeliaz GL37_A_pilon_names2.fasta IRGSP_1.0.genome_names2.fasta
cp sibeliaz_out/* sibeliaz_GL37_A/
cp sibeliaz_GL37_A/alignment.maf GL37_A.maf
ragout GL37_A_recipe -o ragout_GL37_A -s maf --refine --overwrite --repeats --debug -t 30

sibeliaz PDP211_B_pilon_names2.fasta IRGSP_1.0.genome_names2.fasta
cp sibeliaz_out/* sibeliaz_PDP211_B/
cp sibeliaz_PDP211_B/alignment.maf PDP211_B.maf
ragout PDP211_B_recipe -o ragout_PDP211_B -s maf --refine --overwrite --repeats --debug -t 30


sibeliaz SRA2_1_C_pilon_names2.fasta IRGSP_1.0.genome_names2.fasta
cp sibeliaz_out/* sibeliaz_SRA2_1_C/
cp sibeliaz_SRA2_1_C/alignment.maf SRA2_1_C.maf
ragout SRA2_1_C_recipe -o ragout_SRA2_1_C -s maf --refine --overwrite --repeats --debug -t 30


mkdir sibeliaz_J23_A
sibeliaz J23_A_pilon_names2.fasta IRGSP_1.0.genome_names2.fasta
cp sibeliaz_out/* sibeliaz_J23_A/
cp sibeliaz_J23_A/alignment.maf J23_A.maf
ragout J23_A_recipe -o ragout_J23_A -s maf --refine --overwrite --repeats --debug -t 30

###############
#QUALITY CONTROL OF THE 4 GENOMES USING BUSCO FOR COMPLETENESS 
###############

export AUGUSTUS_CONFIG_PATH="/users/sfo503/scratch/chloroplast_round2/busco/Augustus/config/"
busco --augustus -r -i ragout_GL37_A/GL37_A_scaffolds2.fasta -l poales_odb10 -c 20 -o results_poales_GL37_A -m genome
busco --augustus -r -i ragout_PDP211_B/PDP211_B_scaffolds.fasta -l poales_odb10 -c 20 -o results_poales_PDP211_B -m genome
busco --augustus -r -i ragout_SRA2_1_C/SRA2_1_C_scaffolds.fasta -l poales_odb10 -c 20 -o results_poales_SRA2_1_C -m genome
busco --augustus -r -i ragout_J23_A/J23_A_scaffolds.fasta -l poales_odb10 -c 20 -o results_poales_J23_A -m genome

quast.py -t 4 -o quast_GL37_A ragout_GL37_A/GL37_A_scaffolds2.fasta
quast.py -t 4 -o quast_PDP211_B  ragout_PDP211_B/PDP211_B_scaffolds.fasta
quast.py -t 4 -o quast_SRA2_1_C  ragout_SRA2_1_C/SRA2_1_C_scaffolds.fasta
quast.py -t 4 -o quast_J23_A  ragout_J23_A/J23_A_scaffolds.fasta


####################
#COVERAGE of each genome
####################
module load minimap2/2.26-GCCcore-12.2.0 
module load BEDTools/2.30.0-GCC-11.3.0
mkdir coverage
file=names
while read -r line;do
minimap2 -ax map-ont ./ragout/{line}_scaffolds.fasta {line}_normalized1.fq {line}_normalized2.fq  -t SLURM_CPUS_PER_TASK > ./coverage/{line}_mapping_unmapped_raw_to_genome.sam
samtools view -S -bo -h ./coverage/{line}_mapping_unmapped_raw_to_genome.sam -o ./coverage/{line}_mapping_unmapped_raw_to_genome.bam
samtools sort -m 1G ./coverage/{line}_mapping_unmapped_raw_to_genome.bam -o ./coverage/{line}_mapping_unmapped_raw_to_genome_sorted.bam
samtools index ./coverage/{line}_mapping_unmapped_raw_to_genome_sorted.bam
bedtools genomecov -ibam ./coverage/{line}_mapping_unmapped_raw_to_genome_sorted.bam -d > ./coverage/{line}_coverage_genome_base.txt
done < "file"


################
##2. ANNOTATION DE NOVO USING PROTEINS FROM VIRIDIPLANTAE + PROT NIPPONBARE
#DOWNLOAD WITH NIPPONBARE
################
wget -P ./ https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_genome.fasta.gz
wget -P ./ https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_representative_annotation_2023-03-15.tsv.gz
wget -P ./ https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_representative_transcript_exon_2023-03-15.gtf.gz
wget -P ./ https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_representative_2023-03-15.tar.gz
wget -P ./ https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_cds_2023-03-15.fasta.gz
wget -P ./ https://rapdb.dna.affrc.go.jp/download/archive/irgsp1/IRGSP-1.0_protein_2023-03-15.fasta.gz
gzip -d IRGSP*.gz



######BRAKKER for annotations do first RepeatModeler and RepeatMasker for masking and them BRAKER.pl for the de-novo annotation
mkdir ./braker
####masked genome ###
module load BRAKER/2.1.6-foss-2022a
module load RepeatModeler/2.0.4-foss-2022a
module load RepeatMasker/4.1.4-foss-2022a
module load libyaml/0.2.5-GCCcore-11.3.0 
module load PyYAML/6.0-GCCcore-11.3.0 
module load ruamel.yaml/0.17.21-GCCcore-11.3.0
export PATH=PATH:/users/sfo503/perl5
export PERL5LIB=PERL5LIB:/users/sfo503/perl5
export PERL5LIB=PERL5LIB:rice_genomes/lib/perl5/x86_64-linux-thread-multi/YAML
export PATH=libs/ncbi-blast-2.15.0+-src/c++/ReleaseMT/bin:PATH ###instead of this reinstall repeatmodeler and move rmblast to where blast is now which is in libs/


BuildDatabase -name database_GL37_A ./ragout/GL37_A_scaffolds.fasta -engine ncbi
RepeatModeler -database database_GL37_A -engine ncbi -LTRStruct -engine ncbi 
RepeatMasker -gff -lib RM_1549203.TueNov281614322023/consensi.fa ./ragout/GL37_A_scaffolds.fasta -dir GL37_A_masked
BRAKER/scripts/braker.pl --genome GL37_A_masked/GL37_A_scaffolds.fasta.masked --workingdir=braker_GL37_A --hints=rice_genomes/braker/hintsfile.gff --skipAllTraining  --species Sp_2 

########
BuildDatabase -name database_J23_A ./ragout/J23_A_scaffolds.fasta -engine ncbi
RepeatModeler -database database_J23_A -engine ncbi -LTRStruct -engine ncbi 
RepeatMasker -gff -lib rice_genomes/RM_2717714.FriDec10938542023/consensi.fa ./ragout/J23_A_scaffolds.fasta -dir J23_A_masked
BRAKER/scripts/braker.pl --genome J23_A_masked/J23_A_scaffolds.fasta.masked --prot_seq=Viridiplantae_rice.fasta --workingdir=braker_J23

#######
BuildDatabase -name database_SRA2_1_C ./ragout/SRA2_1_C_scaffolds.fasta -engine ncbi
RepeatModeler -database database_SRA2_1_C -engine ncbi -LTRStruct -engine ncbi 
RepeatMasker -gff -lib rice_genomes/RM_88138.FriDec10940342023/consensi.fa ./ragout/SRA2_1_C_scaffolds.fasta -dir SRA2_1_C_masked
BRAKER/scripts/braker.pl --genome SRA2_1_C_masked/SRA2_1_C_scaffolds.fasta.masked --prot_seq=Viridiplantae_rice.fasta --workingdir=braker_SRA2_1_C

######
BuildDatabase -name database_PDP211_B ./ragout/PDP211_B_scaffolds.fasta -engine ncbi
RepeatModeler -database database_PDP211_B -engine ncbi -LTRStruct -engine ncbi 
RepeatMasker -gff -lib rice_genomes/RM_783384.FriDec10941072023/consensi.fa ./ragout/PDP211_B_scaffolds.fasta -dir PDP211_B_masked
BRAKER/scripts/braker.pl --genome PDP211_B_masked/PDP211_B_scaffolds.fasta.masked --prot_seq=Viridiplantae_rice.fasta --workingdir=braker_PDP211_B


###########################
##3.ANNOTATION FILTERING AND FUNCTIONAL CHARACTERIZATION after BRAKER, use EGGNOGG and gFACs
############################
###first eggnog and then filtering with gFACs 
module load BioPerl/1.7.2-GCCcore-8.2.0-Perl-5.28.1 
module load DIAMOND/0.9.30-GCC-8.3.0 

mkdir ./gFACs
mkdir eggnog
mkdir eggnog_database
/users/sfo503/scratch/eggnog-mapper/download_eggnog_data.py -P --data_dir eggnog_database -M  -H -d  33090 -y
export EGGNOG_DATA_DIR=rice_genomes/eggnog_database


Augustus/scripts/gtf2gff.pl <rice_genomes/braker_GL37_A/brakerGL37_A.gtf  --out=rice_genomes/braker_GL37_A/brakerGL37_A.gff 
emapper.py --override --cpu 100 -i rice_genomes/braker_GL37_A/brakerGL37_A.codingseq --output GL37_A --output_dir eggnog -m diamond --itype CDS 

Augustus/scripts/gtf2gff.pl <rice_genomes/braker_J23/brakerJ23_A.gtf  --out=rice_genomes/braker_J23/brakerJ23_A.gff 
emapper.py --override --cpu 100 -i rice_genomes/braker_J23/brakerJ23_A.codingseq --output J23_A --output_dir eggnog -m diamond --itype CDS 

Augustus/scripts/gtf2gff.pl <rice_genomes/braker_PDP211_B/brakerPDP211_B.gtf  --out=rice_genomes/braker_PDP211_B/brakerPDP211_B.gff 
emapper.py --override --cpu 100 -i rice_genomes/braker_PDP211_B/brakerPDP211_B.codingseq --output PDP211_B --output_dir eggnog -m diamond --itype CDS 

Augustus/scripts/gtf2gff.pl <rice_genomes/braker_SRA2_1_C/brakerSRA2_1_C.gtf  --out=rice_genomes/braker_SRA2_1_C/brakerSRA2_1_C.gff 
emapper.py --override --cpu 100 -i rice_genomes/braker_SRA2_1_C/brakerSRA2_1_C.codingseq --output SRA2_1_C --output_dir eggnog -m diamond --itype CDS 


###now can filter with gFAC 
mkdir gFACs
perl /users/sfo503/scratch/gFACs/gFACs.pl -f braker_2.1.5_gtf -p GL37_A_unique --rem-all-incompletes --entap-annotation rice_genomes/eggnog/GL37_A.emapper.annotations --fasta ./ragout/GL37_A_scaffolds.fasta -get-fasta-without-introns --rem-genes-without-start-and-stop-codon --entap-annotation --unique-genes-only --annotated-all-genes-only --create-gtf --get-fasta --distributions exon_lengths intron_lengths CDS_lengths gene_lengths exon_position intron_position --compatibility stringtie_1.3.6_gtf --create-gff3 --get-protein-fasta --statistics -O ./gFACs/ rice_genomes/braker_GL37_A/brakerGL37_A.gtf
perl /users/sfo503/scratch/gFACs/gFACs.pl -f braker_2.1.5_gtf -p J23_A_unique --rem-all-incompletes --entap-annotation rice_genomes/eggnog/J23_A.emapper.annotations --fasta ./ragout/J23_A_scaffolds.fasta -get-fasta-without-introns --rem-genes-without-start-and-stop-codon --entap-annotation --unique-genes-only --annotated-all-genes-only --create-gtf --get-fasta --distributions exon_lengths intron_lengths CDS_lengths gene_lengths exon_position intron_position --compatibility stringtie_1.3.6_gtf --create-gff3 --get-protein-fasta --statistics -O ./gFACs/ rice_genomes/braker_J23/brakerJ23_A.gtf
perl /users/sfo503/scratch/gFACs/gFACs.pl -f braker_2.1.5_gtf -p PDP211_B_unique --rem-all-incompletes --entap-annotation rice_genomes/eggnog/PDP211_B.emapper.annotations --fasta ./ragout/PDP211_B_scaffolds.fasta -get-fasta-without-introns --rem-genes-without-start-and-stop-codon --entap-annotation --unique-genes-only --annotated-all-genes-only --create-gtf --get-fasta --distributions exon_lengths intron_lengths CDS_lengths gene_lengths exon_position intron_position --compatibility stringtie_1.3.6_gtf --create-gff3 --get-protein-fasta --statistics -O ./gFACs/ rice_genomes/braker_PDP211_B/brakerPDP211_B.gtf
perl /users/sfo503/scratch/gFACs/gFACs.pl -f braker_2.1.5_gtf -p SRA2_1_C_unique --rem-all-incompletes --entap-annotation rice_genomes/eggnog/SRA2_1_C.emapper.annotations --fasta ./ragout/SRA2_1_C_scaffolds.fasta -get-fasta-without-introns --rem-genes-without-start-and-stop-codon --entap-annotation --unique-genes-only --annotated-all-genes-only --create-gtf --get-fasta --distributions exon_lengths intron_lengths CDS_lengths gene_lengths exon_position intron_position --compatibility stringtie_1.3.6_gtf --create-gff3 --get-protein-fasta --statistics -O ./gFACs/ rice_genomes/braker_SRA2_1_C/brakerSRA2_1_C.gtf

####DO NOW FINAL EGGNOG AFTER gFACs
mkdir eggnog_after_gFACs
emapper.py --override --cpu 100 -i ./gFACs/GL37_A_unique_genes.fasta --output GL37_A_after_gFACs --output_dir eggnog_after_gFACs -m diamond --itype CDS --decorate_gff yes
emapper.py --override --cpu 100 -i ./gFACs/J23_A_unique_genes.fasta --output J23_A_after_gFACs --output_dir eggnog_after_gFACs -m diamond --itype CDS --decorate_gff yes
emapper.py --override --cpu 100 -i ./gFACs/PDP211_B_unique_genes.fasta --output PDP211_B_after_gFACs --output_dir eggnog_after_gFACs -m diamond --itype CDS --decorate_gff yes
emapper.py --override --cpu 100 -i ./gFACs/SRA2_1_C_unique_genes.fasta --output SRA2_1_C_after_gFACs --output_dir eggnog_after_gFACs -m diamond --itype CDS --decorate_gff yes




##################
##4.EDTA to annotate transposable elements
#################


#####EDTa
EDTA.pl --genome ./ragout/GL37_A_scaffolds.fasta --species Rice --step all --cds rice_genomes/braker_GL37_A/brakerGL37_A.codingseq --sensitive 1 --anno 1 --evaluate 0 --threads 40 --overwrite 1 --force 1 
EDTA.pl --genome ./ragout/J23_A_scaffolds.fasta --species Rice --step all --cds rice_genomes/braker_J23/brakerJ23_A.codingseq --sensitive 1 --anno 1 --evaluate 0 --threads 40 --overwrite 1 --force 1 
EDTA/EDTA.pl --genome ./ragout/PDP211_B_scaffolds.fasta --species Rice --step all --cds rice_genomes/braker_PDP211_B/brakerPDP211_B.codingseq --sensitive 1 --anno 1 --evaluate 0 --threads 40 --overwrite 1 --force 1 
EDTA/EDTA.pl --genome ./ragout/SRA2_1_C_scaffolds.fasta --species Rice --step all --cds rice_genomes/braker_SRA2_1_C/brakerSRA2_1_C.codingseq --sensitive 1 --anno 1 --evaluate 0 --threads 40 --overwrite 1 --force 1 
####species can be others





####################
##5.MASH DISTANCE to find how different the 4 genomes are between them. MASH gives dissimilarity 
####################
module load  Mash/2.2-foss-2019b
mash sketch rice_genomes/ragout/GL37_A_scaffolds.fasta
mash sketch rice_genomes/ragout/J23_A_scaffolds.fasta
mash sketch rice_genomes/ragout/PDP211_B_scaffolds.fasta
mash sketch rice_genomes/ragout/SRA2_1_C_scaffolds.fasta

mash paste mash_rice -l mash_rice 
mash triangle -l mash_rice -E -k 11 > mash_results
mash triangle -l mash_rice -E > mash_results2




####################
##6. PHYLOGENY USING THE GENE PHYTOCHROME C (Os03g0752100-01) AND DOWNLOADING SEQUENCES FROM https://agrigenome.dna.affrc.go.jp/tasuke/ricegenomes/
####################
#####get the sequence corresponding with the PhyC using as input 
file=rice_names ###file with the 4 genomes names GL37_A, J23_A, PDP211_B, SRA2_1_C
while read -r line;do
makeblastdb -in {line}_scaffolds.fasta -out {line}_scaffolds.fasta -dbtype nucl -title {line}_scaffolds -parse_seqids
blastn -query phylogeny/phyC -db /mnt/scratch/users/sfo503/rice_genomes/{line}_scaffolds.fasta -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq"  -max_target_seqs 1 -max_hsps 1 > ./phylogeny/{line}_phyC
awk '{print ">" 2,9,10"\n"13}'  ./phylogeny/{line}_phyC  > ./phylogeny/{line}_phyC.fa ####column 9 and 10 are the coordinates and column 2 is the contig
sed 's/>/>_/g' ./phylogeny/{line}_phyC.fa > ./phylogeny/{line}_phyC_0.fa
sed 's/ /_/g' ./phylogeny/{line}_phyC_0.fa > ./phylogeny/{line}_phyC_1.fa
done < "file"

mv *_1.fa ./MSA
####rename the fasta sequences with the name in the file
for file in *.fa;
   do
       sed -i "s/>/>{file%%.*}/" "file" ;
done


module load SeqKit/2.3.1
##remove dups ####
seqkit rmdup -s < TASUKE.fa > TASUKE2.fa
sed 's/ /_/g' TASUKE2.fa > TASUKE_out.fa
#####this is to add a number  to the sequences from the database that are repeated
awk 'BEGIN{RS=">";OFS="\n"}(NR>1){print ">"1"_"(NR-1)"\n";1="";print 0}' TASUKE_out.fa | awk '0' > ./MSA/TASUKE_out.fa

###generate tree
cat ./*.fa > fasta_combined.fasta
mafft --thread 40 fasta_combined.fasta > fasta_combined_aligned
trimal -in fasta_combined_aligned -out fasta_combined_aligned_trimmed_output_fasta -automated1 -fasta ####checked manually and removed part. the total length now is 2357bp
iqtree -s fasta_combined_aligned_trimmed_output_manual.fas --alrt 1000 -B 1000

