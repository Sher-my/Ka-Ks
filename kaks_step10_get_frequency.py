#!/bin/python
from Bio import SeqIO
import numpy as np
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
from statsmodels.stats.multicomp import pairwise_tukeyhsd
path = "C:/summary/kaks/motif/"
file_names = ['4dtv_homo_seq', '4dtv_upstream_seq', 'GRCh38_random', 'GRCh38_random1', 'ka_homo_seq', 'ka_upstream_seq']
for fn in file_names:
        # 读取fasta文件中的所有序列，并获取第一条序列以确定公共长度
        fasta_sequences = SeqIO.parse(path + fn + ".fa", "fasta")
        first_sequence = next(fasta_sequences)
        common_length = len(first_sequence)

        # 初始化一个numpy数组用于存储各位置的碱基频次
        base_counts = np.zeros((common_length, 4), dtype=int)
        base_labels = ["A", "C", "G", "T"]

        # 遍历所有序列并统计碱基频率
        for sequence in fasta_sequences:
            if len(sequence) != common_length: 
                print("Error: Not all sequences are of the same length!")
                break
            else:
                # 将DNA序列转换为字符串以方便计数
                seq_str = str(sequence.seq).upper()
                
                for pos, base in enumerate(seq_str):
                    base_index = base_labels.index(base)
                    base_counts[pos, base_index] += 1

        # 创建DataFrame形式的结果，并添加"position"列
        freq_table = pd.DataFrame(base_counts, columns=base_labels)
        freq_table.insert(0, "position", range(1, common_length + 1))
        # freq_table.to_csv(path + fn + ".csv", index=False)
        # # 提取碱基频次矩阵
        freq_matrix = freq_table.iloc[:, 1:].values
        # print(freq_matrix)
        # # 创建Logo对象
        # logo_df = pd.DataFrame(freq_matrix, columns=base_labels, index=freq_table['position'])
        # # 假设您已经有了一个DataFrame（logo_df）, 其中列是碱基，行是位置，值是每个位置上碱基的计数或频率
        # # 创建Logo对象并设置字体高度范围
        # logo = logomaker.Logo(logo_df)  # 设置stack_width以控制碱基间的间距
        # logo.height_scale = 'auto'  # 这里使用'auto'选项，它会根据数据自适应地调整logo的高度

        # # 调整整个logo的全局缩放比例
        # logo.figure_size = (8, 1)  # 可以调整figure_size参数来改变整体大小

        # # 绘制logo
        # logo.draw()

        # # 显示或保存图像
        # plt.show()

        # 将各列数据转换为系列列表
        data = [freq_matrix[:, i] for i in range(4)]

        # 绘制箱型图
        import matplotlib.pyplot as plt

        plt.boxplot(data, labels=['A', 'C', 'G', 'T'])
        plt.xlabel('Nucleotide')
        plt.ylabel('Frequency')
        plt.title('Nucleotide Frequency Distribution')
        plt.savefig(path + fn + '_boxplot.pdf')
        plt.clf()