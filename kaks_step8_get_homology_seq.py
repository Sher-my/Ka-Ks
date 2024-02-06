#!/bin/python
import pandas as pd
def extract_sequences(file1_path, file2_path, output_path):
    df1 = pd.read_csv(file1_path, sep='\t', header=0)
    fasta_dict = {}
    with open(file2_path) as f:
        sequence = ''
        for line in f:
            if line.startswith('>'):
                if sequence:
                    fasta_dict[person_name] = sequence.strip()
                person_name = line[1:].strip()
                sequence = ''
            else:
                sequence += line.strip()

        fasta_dict[person_name] = sequence # 添加最后一个序列到字典中（防止遗漏）

    with open(output_path, 'w') as out_f: # 创建一个新的fasta文件，并写入提取的子序列
        for _, row in df1.iterrows():
            person_name = row['Person_Name']
            start = int(row['start'])
            end = int(row['end'])
            
            if person_name in fasta_dict: # 确保ID存在于fasta字典中
                seq = fasta_dict[person_name][start:end]
                out_f.write(f'>{person_name}\n{seq}\n')
            else:
                print(f"Warning: '{person_name}' not found in the FASTA file.")

files = ['filtered_4dtv_corrected_data.tsv', 'filtered_ka_ks_data.tsv']
for file in files:  
    extract_sequences(file, 'GRCh38_one.cds', file.split('_')[1] + '_homo_seq.fa')