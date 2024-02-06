#!/bin/python
from Bio import SeqIO
import random

def get_valid_subseq(seq, length=60):
    # 确保序列长度至少是length的3倍，因为要保证子序列能完整翻译
    if len(seq) < 3 * length:
        return None
    
    # 随机选取起始点，保证子序列长度是3的倍数
    start = random.randint(0, len(seq) - 3 * length)
    
    # 提取满足条件的子序列
    sub_seq = seq[start : start + 3 * length]
    
    return sub_seq

# 读取原始fasta文件
with open("GRCh38_one.cds", "r") as handle:
    records = list(SeqIO.parse(handle, "fasta"))

# 创建一个新的fasta记录列表用于存储随机子序列
new_records = []

# 对于每条序列
for record in records:
    sub_seq = get_valid_subseq(record.seq, 20)  # 这里以20个密码子对应长度60的子序列为例

    # 如果成功获取到符合条件的子序列，则添加至新记录列表
    if sub_seq is not None:
        new_record = SeqIO.SeqRecord(sub_seq, id=f"{record.id}_random_subseq", description="")
        new_records.append(new_record)

# 将新的fasta记录写入输出文件
with open("GRCh38_random.fa", "w") as output_handle:
    SeqIO.write(new_records, output_handle, "fasta")