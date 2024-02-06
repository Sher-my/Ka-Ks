#!/bin/python
import argparse
parser = argparse.ArgumentParser()
parser.description='Please enter specie and length of aa...'
parser.add_argument("-s", "--specie", dest='s', help="enter a str like 'GRCh38'", type=str, default='GRCh38')
parser.add_argument("-l", "--length", dest='l', help="enter an int like '10'", type=int, default=10)
args = parser.parse_args()
specie = args.s
length_aa = args.l
def remove_duplicate_sequences(file_in, file_out):
    sequence_dict = {}

    with open(file_in, 'r') as f:
        for line in f:
            if line.startswith(">"):
                seq_name = line[1:].strip()
                next_line = next(f).strip()
                
                if next_line in sequence_dict: # 如果序列已存在于字典中，则添加到名称列表；否则新建一个名称列表
                    sequence_dict[next_line].append(seq_name)
                else:
                    sequence_dict[next_line] = [seq_name]

    with open(file_out, 'w') as f: # 写入新文件，合并序列名并保留不重复的序列
        for seq, names in sequence_dict.items():
            combined_name = "&".join(names) # 使用&符号连接序列名
            f.write(f">{combined_name}\n{seq}\n")

remove_duplicate_sequences('/disk/yt/kaks/cds/' + specie + '_' + length_aa + 'aa.cds', '/disk/yt/kaks/cds/' + specie + '_' + length_aa + 'aa_rmdup.cds')