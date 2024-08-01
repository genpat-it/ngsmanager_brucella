#!/usr/bin/env python3

import argparse
import csv
import re
import subprocess


def prepare_input(file_list, file_aliases, reference, output):
    aliases = {}
    for line in csv.reader(file_aliases, delimiter="\t"):
        if line[2] == 'SILAB2':
            aliases[line[0]] = line[1]
    subprocess.run("awk 1 " + reference + " > " + output, shell=True, check=True)
    for line in csv.reader(file_list, delimiter="\t"):
        filename = line[0]
        match = re.search('DS\d+-DT\d+_([^_]+)_.*\.fa\w*', filename)
        if not match:
            raise Exception("could not extract sample code from file name: " + filename)
        sample = match.group(1)
        alias = aliases.get(sample, sample)
        anno = alias.replace("FROM-", "").split(".")[0]
        if alias[-4:] != ".1.1":
            fromm = alias.split(".")[2] + "." + alias.split(".")[3] + "." + alias.split(".")[4]
        else:
            fromm = alias.split(".")[2]
        with open(filename, 'r') as current_fasta:
            header = current_fasta.readline().strip()
            cmd = "sed 's~{}~>hCoV-19/Italy/ABR-IZSGC-{}/{}~' {} > consensus.tmp".format(header, fromm, anno, filename)
            subprocess.run(cmd, shell=True, check=True)
            subprocess.run("awk 1 consensus.tmp >> {}".format(output), shell=True, check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--file_list', dest='file_list', required=True,
                        help='file with the list of consensus to be processed',
                        type=argparse.FileType('r'))
    parser.add_argument('--aliases', dest='file_aliases', required=True,
                        help='aliases file',
                        type=argparse.FileType('r'))
    parser.add_argument('--reference', dest='reference', required=True,
                        help='reference file',
                        type=str)
    parser.add_argument('--output', dest='output', required=True,
                        help='multi-fasta file name',
                        type=str)
    args = parser.parse_args()

    prepare_input(args.file_list, args.file_aliases, args.reference, args.output)
