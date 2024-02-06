#!/bin/python
import csv

input_file = 'all-results.txt'

ka_kas_output = 'ka_ks_data.tsv'
fourDTV_output = '4dtv_corrected_data.tsv'
filtered_ka_kas_output = 'filtered_ka_ks_data.tsv'
filtered_fourDTV_output = 'filtered_4dtv_corrected_data.tsv'

with open(input_file, 'r') as f:
    reader = csv.reader(f, delimiter=' ')
    header = next(reader)  # 跳过标题行
    seq_index = header.index('Seq')
    ka_ks_index = header.index('Ka_Ks')
    fourDTV_corr_index = header.index('4dtv_corrected')

    data = []
    for row in reader:
        mouse_humans = row[seq_index].split('-')[1].split('&')
        ka_ks_val = row[ka_ks_index]
        fourDTV_corr_val = row[fourDTV_corr_index]

        # 处理'NA'值
        if ka_ks_val == 'NA':
            ka_ks_val = None
        else:
            ka_ks_val = float(ka_ks_val)

        if fourDTV_corr_val == 'NA':
            fourDTV_corr_val = None
        else:
            fourDTV_corr_val = float(fourDTV_corr_val)

        for mouse_human in mouse_humans:
            try:
                person_name, person_id = mouse_human.split('+')
                ka_kas_line = [person_name, person_id, int(person_id) + 60, ka_ks_val]
                fourDTV_line = [person_name, person_id, int(person_id) + 60, fourDTV_corr_val]

                # 只包含非'NA'记录
                if ka_ks_val is not None and fourDTV_corr_val is not None:
                    data.append((person_name, ka_kas_line, fourDTV_line))
            except ValueError:
                continue  # 如果分割得到的不是两个部分，则跳过此对应关系

# 写入原始处理后的数据（只包含非'NA'记录）
with open(ka_kas_output, 'w', newline='') as kf, open(fourDTV_output, 'w', newline='') as ff:
    ka_kas_writer = csv.writer(kf, delimiter='\t')
    fourDTV_writer = csv.writer(ff, delimiter='\t')

    ka_kas_writer.writerow(['Person_Name', 'start', 'end', 'Ka_Ks'])
    fourDTV_writer.writerow(['Person_Name', 'start', 'end', '4dtv_corrected'])

    for name, ka_kas_info, fourDTV_info in data:
        ka_kas_writer.writerow(ka_kas_info)
        fourDTV_writer.writerow(fourDTV_info)

# 筛选显著值（计算每个名称下的平均值并选择大于平均值的记录）
significant_data = {}
for person, ka_kas_info, fourDTV_info in data:
    if person not in significant_data:
        significant_data[person] = {
            'ka_ks': {'all_values': [], 'max_value': None, 'lines': []},
            'fourDTV': {'all_values': [], 'max_value': None, 'lines': []}
        }
    if ka_kas_info[3] is not None:
        significant_data[person]['ka_ks']['all_values'].append(ka_kas_info[3])
    if fourDTV_info[3] is not None:
        significant_data[person]['fourDTV']['all_values'].append(fourDTV_info[3])

# 计算平均值并筛选数据
for person, info in significant_data.items():
    ka_ks_avg = sum(info['ka_ks']['all_values']) / len(info['ka_ks']['all_values']) if info['ka_ks']['all_values'] else None
    for ka_kas_line in (row[1] for row in data if row[0] == person and row[1][3] is not None and row[1][3] > ka_ks_avg):
        info['ka_ks']['lines'].append(ka_kas_line)

    fourDTV_avg = sum(info['fourDTV']['all_values']) / len(info['fourDTV']['all_values']) if info['fourDTV']['all_values'] else None
    for fourDTV_line in (row[2] for row in data if row[0] == person and row[2][3] is not None and row[2][3] < fourDTV_avg):
        info['fourDTV']['lines'].append(fourDTV_line)

# 写入筛选后的显著值数据（这里假设每个名称至少有一个有效值大于平均值，实际情况需要处理空值情况）
with open(filtered_ka_kas_output, 'w', newline='') as fkf, open(filtered_fourDTV_output, 'w', newline='') as fff:
    filtered_ka_kas_writer = csv.writer(fkf, delimiter='\t')
    filtered_fourDTV_writer = csv.writer(fff, delimiter='\t')

    filtered_ka_kas_writer.writerow(['Person_Name', 'start', 'end', 'Ka_Ks'])
    filtered_fourDTV_writer.writerow(['Person_Name', 'start', 'end', '4dtv_corrected'])

    for person, info in significant_data.items():
        for line in info['ka_ks']['lines']:
            filtered_ka_kas_writer.writerow(line[:3] + [line[3]])
        for line in info['fourDTV']['lines']:
            filtered_fourDTV_writer.writerow(line[:3] + [line[3]])