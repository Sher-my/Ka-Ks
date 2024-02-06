#!/bin/python
import csv
from Bio import SeqIO

ts = ['MS', 'ribo', 'statistical']
for t in ts:
    file1_path = 'pep_hg38.txt'
    fasta_in_path = t + '_CR.fa'
    fasta_out_path = t + '_CR_add.fa'

    motif_data = {}
    with open(file1_path, 'r') as f1:
        reader = csv.reader(f1, delimiter='\t')
        header = next(reader)  # 跳过标题行
        for row in reader:
            id_, motif, start, stop = row[0], row[1], row[2], row[3]
            motif_data[id_] = f"{motif}+{start}+{stop}"

    with open(fasta_in_path, 'r') as fi, open(fasta_out_path, 'w') as fo:
        for record in SeqIO.parse(fi, 'fasta'):
            id_ = record.id
            if id_ in motif_data:
                new_id = f"{id_}+{motif_data[id_]}"
                record.id = new_id
                record.description = ""  # 清空描述以避免与序列名中的信息冲突
                SeqIO.write(record, fo, 'fasta')
            else:
                print(f"Warning: ID {id_} not found in file1.")
                SeqIO.write(record, fo, 'fasta')