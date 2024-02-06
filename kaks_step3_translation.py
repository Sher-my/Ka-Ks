#!/bin/python
import argparse
import re

def translate_cds_to_pep(cds_path, pep_path):
    path = "/disk/yt/kaks/"
    parser = argparse.ArgumentParser(description='Please enter input and output file names...')
    parser.add_argument("-fni", "--file_name_in", dest='fni', help="enter a string like 'GRCh38_Naa_rmdup.cds'", type=str, default='GRCh38_20aa_rmdup.cds')
    parser.add_argument("-fno", "--file_name_out", dest='fno', help="enter a string like 'GRCh38_Naa.pep'", type=str, default='GRCh38_20aa.pep')
    args = parser.parse_args()
    specie_fi = args.fni
    specie_fo = args.fno

    codon_table = {
            'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
            'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
            'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
            'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
            'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
            'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
            'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
            'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
            'UAC': 'Y', 'UAU': 'Y', 'UAA': '*', 'UAG': '*',
            'UGC': 'C', 'UGU': 'C', 'UGA': '*', 'UGG': 'W',
            }
    
    cds_file = f"{path}cds/{specie_fi}"
    pep_file = f"{path}protein/{specie_fo}"

    with open(cds_file, 'r') as cds, open(pep_file, 'w') as out:
        cds_records = [line.strip() for line in cds.readlines()]
        for i in range(0, len(cds_records), 2):
            id_line = cds_records[i]
            rna_seq = cds_records[i + 1].replace('T', 'U')

            protein_seq = ''.join([codon_table[codon] for codon in [rna_seq[j:j+3] for j in range(0, len(rna_seq), 3) if len(rna_seq[j:j+3]) == 3]])

            out.write(id_line)
            out.write(protein_seq + '\n')

if __name__ == "__main__":
    translate_cds_to_pep()