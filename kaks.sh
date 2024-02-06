#!/bin/bash
###WGD analysis(ka/ks and 4Dtv)
###prepare files:specie.cds.fa/specie.protein.fa/specie.gff3

# Set variables
SPECIES=("GRCm39" "GRCh38")
N=20
OUTPUT_DIR="01blast"

# Prepare one-line cds files
for specie in "${SPECIES[@]}"; do
  awk '/>/ {name=$0; seq[name]=seq[name]RS$0} END{for(i in seq) print i ORS seq[i]}' "cds/${specie}.cds" >> "cds/${specie}_one.cds"
done

# Run kaks pipeline steps
python kaks_step1_get_refrence_seq_Naa.py
python kaks_step2_rmdup_seq.py
python kaks_step3_translation.py

# Activate blast environment and create blast directory
#conda activate blast
mkdir -p "$OUTPUT_DIR"

# Step 1: Blastp
makeblastdb -in "protein/GRCm39_${N}aa.pep" -dbtype prot -out "${OUTPUT_DIR}/GRCm39_${N}aa"
blastp -query "protein/GRCh38_${N}aa.pep" -out "${OUTPUT_DIR}/GRCh38_GRCm39.blast" -db "${OUTPUT_DIR}/GRCm39_${N}aa" -outfmt 6 -evalue 1e-5 -num_alignments 5 > "${OUTPUT_DIR}/blastp.log" 2>&1

# Step 2: GFF extraction (method 1)
for specie in "${SPECIES[@]}"; do
  blastn -db "fa/index/${specie}" -query "cds/${specie}_${N}aa_rmdup.cds" -out "fa/blast/${specie}.txt" -outfmt '6 pident sseqid qseqid sstart send' -max_hsps 1 -num_alignments 1
done
#conda activate base
python fa/blast/kaks_step4_get_gff.py
mv fa/blast/*.gff "$OUTPUT_DIR"

# Step 2 (alternative): Convert gff3 to gff (method 2)
awk '$3=="transcript"{print $1"\t"$4"\t"$5"\t"$9}' GRCh38.gff3 | awk -F 'ID=' '{print $1$2}' | awk -F 'Parent=' '{print $1}' | awk -F ';' '{print $1"\t"$4"\t"$2"\t"$3}' > GRCh38.gff

# Step 3: MCScanX
cd "$OUTPUT_DIR"
/disk/yt/biotools/MCScanX/MCScanX GRCh38_GRCm39 -g -3 -e 1e-10
grep "ENS" GRCh38_GRCm39.collinearity | awk '{print $3"\t"$4}' > GRCh38_GRCm39.homolog

# Step 4: Ka/Ks and 4Dtv calculation
SPECIE="GRCh38"
GROUP="GRCm39"
echo "12" > "${SPECIE}_${GROUP}.proc"
/disk/yt/biotools/ParaAT2.0/ParaAT.pl -h "${OUTPUT_DIR}/GRCh38_GRCm39.homolog" -n "cds/${SPECIE}_${GROUP}.cds" -a "protein/${SPECIE}_${GROUP}.pep" -m clustalw2 -p "${SPECIE}_${GROUP}.proc" -f axt -o "${SPECIE}_${GROUP}_out" 2> ParaAT.log

cd "${SPECIE}_${GROUP}_out"
ls > list.txt

while read i; do
    KaKs_Calculator -i "$i" -o "${i}.kaks" -m YN
    python ../axt2one-line.py "$i" "${i}.one-line"
    perl ../calculate_4DTV_correction.pl "${i}.one-line" >"${i}.4dtv"
    awk 'NR>1{print $1"\t"$3}' "${i}.4dtv" >>all-4dtv.txt
    awk 'NR>1{print $1"\t"$3"\t"$4"\t"$5}' "${i}.kaks" >>all-kaks.txt
done <list.txt

sort all-4dtv.txt | uniq >all-4dtv.results
sort all-kaks.txt | uniq >all-kaks.results
join -1 1 -2 1 all-kaks.results all-4dtv.results >all-results.txt
sed -i '1i(Seq\tKa\tKs\tKa/Ks\t4dtv_corrected)' all-results.txt

cd ..
mv "${SPECIE}_${GROUP}_out"/all-results.txt ./

### step5 get significant_value from kaks&4dtv
python kaks_step5_get_siginificant_value.py

### step6 get CR position
python kaks_step6_get_enid_from_CR.py
#conda activate blast
makeblastdb -in ../cds/GRCh38_one.cds -dbtype nucl -parse_seqids -out GRCh38
blastn -db GRCh38 -query statistical_CR_add.fa -out statistical_CR_add.txt -outfmt ' 6 pident sseqid qseqid sstart send ' -max_hsps 1  -num_alignments 1
blastn -db GRCh38 -query ribo_CR_add.fa -out ribo_CR_add.txt -outfmt ' 6 pident sseqid qseqid sstart send ' -max_hsps 1  -num_alignments 1
blastn -db GRCh38 -query MS_CR_add.fa -out MS_CR_add.txt -outfmt ' 6 pident sseqid qseqid sstart send ' -max_hsps 1  -num_alignments 1

### step7 count kaks/4dtv of CR
python kaks_step7_count_result1.py

### step8 get_seq
python kaks_step8_get_homology_seq.py
python kaks_step9_get_updown_seq.py

### step9 get random_seq
python kaks_step9.1_get_random_seq.py

### step10 get motif
python kaks_step10_get_frequency.py
