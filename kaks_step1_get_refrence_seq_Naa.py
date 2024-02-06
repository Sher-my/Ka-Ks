#!/bin/python
import argparse

def main():
    parser = argparse.ArgumentParser(description='Please enter species and length of amino acids...')
    parser.add_argument("-s", "--species", dest='species', help="enter a string like 'GRCh38'", type=str, default='GRCh38')
    parser.add_argument("-l", "--length_aa", dest='length_aa', help="enter an integer like '10' for aa length", type=int, default=10)
    
    args = parser.parse_args()
    species = args.species
    length_aa = args.length_aa
    
    cds_file_path = f"/disk/yt/kaks/cds/{species}_one.cds"
    output_file_path = f"/disk/yt/kaks/cds/{species}_{length_aa}aa.cds"

    cds_dict = {}
    new_dict = {}

    with open(cds_file_path, 'r') as cds:
        for line in cds:
            if line.startswith('>'):
                id = line.strip()[1:]
                seq = next(cds).strip()
                cds_dict[id] = seq

    for key, seq in cds_dict.items():
        for i in range(0, len(seq), length_aa * 3):
            new_seq = seq[i:i + length_aa * 3]
            if len(new_seq) >= length_aa * 3:
                new_dict[f"{key}+{i}"] = new_seq

    with open(output_file_path, 'w') as out:
        for key, value in new_dict.items():
            out.write(f">{key}\n{value}\n")

if __name__ == "__main__":
    main()