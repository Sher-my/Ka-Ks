#!/bin/python
import pandas as pd
import numpy as np

def count(file_name, file2, col_name):
    df1 = pd.read_table(file_name, sep='\t', header=0)

    file = file2
    df2 = pd.read_table(file, sep='\t', header=None, names= ['Person_Name', 'CR_start', 'CR_end'])

    with open('GRCh38_one.cds') as f:
        seq_dict = {}
        id = None
        sequence = ''
        for line in f:
            if line.startswith('>'):
                if id:
                    seq_dict[id] = sequence
                id = line[1:].strip()
                sequence = ''
            else:
                sequence += line.strip()

    # 将序列长度加入到df1和df2中
    df1['seq_length'] = df1['Person_Name'].map(lambda x: len(seq_dict[x]))
    df2['seq_length'] = df2['Person_Name'].map(lambda x: len(seq_dict[x]))

    # 去重
    df1 = df1.drop_duplicates()
    df2 = df2.drop_duplicates()

    merged_df = pd.merge(df1, df2, on='Person_Name')

    def is_in_range(row):
        return row['start'] >= row['CR_end'] and row['end'] <= row['CR_end'] + 150

    def random_sequence_overlap(row):
        seq_len = row['seq_length_x']
        random_start = np.random.randint(0, seq_len - 150)
        return row['start'] >= random_start and row['end'] <= random_start + 150

    merged_df['in_FS'] = merged_df.apply(is_in_range, axis=1)
    merged_df['random_range'] = merged_df.apply(random_sequence_overlap, axis=1)
    merged_df['out_FS'] = merged_df['random_range'] * merged_df[col_name]

    # 保留要么在FS里为Ture要么在random里为Ture
    filtered_df = merged_df[(merged_df['in_FS'] | merged_df['random_range'])]

    filtered_df[['in_FS', col_name]].to_csv(file_name.split('.')[0] + '&' + file2.split('.')[0] + '&' + col_name.split('.')[0] + '_output.tsv', sep='\t', index=False)
    filtered_df[['random_range', col_name]].to_csv(file_name.split('.')[0] + '&' + file2.split('.')[0] + '&' + col_name.split('.')[0] + '_output_random.tsv', sep='\t', index=False)

    import matplotlib.pyplot as plt
    from scipy.stats import ttest_ind

    # 提取两组数据
    in_fs_values = filtered_df[col_name][filtered_df['in_FS']== True]
    out_fs_values = filtered_df[col_name][filtered_df['random_range']== True]
    
    # 绘制箱型图
    plt.figure(figsize=(10, 6))
    plt.boxplot([in_fs_values, out_fs_values], labels=['in_FS', 'random_range'])
    plt.ylabel(col_name + 'values')
    plt.title('Distribution of ' + col_name + ' Values in Both Conditions')
    # plt.axhline(y=in_fs_values.median(), color='r', linestyle='--', label='Median (in_FS)')
    # plt.axhline(y=out_fs_values.median(), color='g', linestyle='--', label='Median (random_range)')
    plt.legend()

    # 进行显著性检验
    t_statistic, p_value = ttest_ind(in_fs_values, out_fs_values, equal_var=False)

    # 在图形上标注p值
    plt.annotate(f'p-value: {p_value:.3f}', xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top')

    # 保存图片
    plt.savefig(file_name.split('.')[0] + '&' + file2.split('.')[0] + '&' + col_name.split('.')[0] + '_boxplot.pdf')
    plt.show()

filenames = ['filtered_4dtv_corrected_data.tsv', 'filtered_ka_ks_data.tsv']
file2s = ['MS_CR_position.txt', 'ribo_CR_position.txt', 'statistical_CR_position.txt']
for fn in filenames:
    for f2 in file2s:
        if '4dtv' in fn:
            cn = '4dtv_corrected'
        if 'ka' in fn:
            cn = 'Ka_Ks'
        count(fn, f2, cn)