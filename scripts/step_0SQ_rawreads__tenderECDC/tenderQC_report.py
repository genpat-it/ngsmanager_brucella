#!/usr/bin/env python3

import argparse
import multiprocessing as mp
import subprocess

from Bio import SeqIO
from numpy import average

genomesizes = {
               'Acinetobacter': 4,
               'Campylobacter': 2,
               'Clostridium': 4,
               'ENTEROBACTERIACEAE': 5,
               'Escherichia': 5,
               'Legionella' : 3.5,
               'Listeria': 3, 
               'Neisseria': 2,
               'Salmonella': 5,
               'Tuberculosis': 4.5,
               "coronavirus": 0.03
            }


def safe_open(file, mode='rt'):
    # safe open file
    if file.endswith('.gz'):
        import gzip
        return gzip.open(file, mode)
    else:
        return open(file, mode)


def reads_info(fastq):
    res = {}
    len_reads = []
    reads_qual30 = 0
    avg_qual = []
    handle = safe_open(fastq)
    for record in SeqIO.parse(handle, "fastq"):
        len_reads.append(len(record.seq))
        phred_quality_list = record.letter_annotations['phred_quality']
        avg_qual.append(average(phred_quality_list))
        if average(phred_quality_list) > 30:
            reads_qual30 += 1
    handle.close()

    q30_reads = round(reads_qual30 / float(len(len_reads)) * 100, 2)
    mbases = round(float(sum(len_reads)) / 1000000, 2)

    res.update({"avgLengthReads": round(average(len_reads), 2)})
    res.update({"avg_qual": avg_qual})
    # res.update({"len_reads":len_reads})
    # res.update({"reads_qual30":reads_qual30})
    res.update({"q30_reads": q30_reads})
    res.update({"mbases": mbases})
    res.update({"totReads": len(len_reads)})
    return res


def check_contamination(genus_data, sample_info):
    g = open(genus_data)
    check = ""
    for genus_data_row in g.readlines():
        name, taxonomy_id, taxonomy_lvl, kraken_assigned_reads, added_reads, new_est_reads, fraction_total_reads = \
            genus_data_row.strip().split("\t")
        if name in sample_info:
            if float(fraction_total_reads) < 0.95:
                check = "FAIL"
            else:
                check = "PASS"
            break
    if check == "":
        check = "NA"
    return check, int(new_est_reads), float(fraction_total_reads)


def calc_cov(sample_info,assigned,r1_avgLengthReads,r2_avgLengthReads):
    cov = 0
    check = ""
    for s in genomesizes:
        if s in sample_info:
            g_size = genomesizes[s]
            assigned_bases = (assigned*r1_avgLengthReads)+(assigned*r2_avgLengthReads)
            cov = assigned_bases / (g_size*1000000)
            if cov < 50:
                check = "FAIL"
            else:
                check = "PASS"
            break
    return cov, check


def calc_md5(file):
    cmd = "md5sum " + file + " | awk '{print $1}' "
    md5 = subprocess.check_output(cmd, shell=True)
    return md5.strip()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='FASTQ Report')
    parser.add_argument('-n', '--name', help='Sample name', type=str, required=True)
    parser.add_argument('-r1', '--reads1', help='Fastq R1', type=str, required=True)
    parser.add_argument('-r2', '--reads2', help='Fastq R2', type=str, required=True)
    parser.add_argument('-g', '--genus', help='Kraken genus', type=str, required=True)
    parser.add_argument('-sp', '--species', help='Species to check', type=str, required=True)
    args = parser.parse_args()

    pool = mp.Pool(processes=2)
    fastqInfo = pool.map(reads_info, [args.reads1, args.reads2])
    r1 = fastqInfo[0]
    r2 = fastqInfo[1]
    sampleName = args.name
    sampleSpecies = args.species
    checkGenus,assigned, fraction = check_contamination(args.genus, sampleSpecies)
    coverage, checkCov = calc_cov(sampleSpecies, assigned,r1["avgLengthReads"],r2["avgLengthReads"])
    r1_md5 = calc_md5(args.reads1).decode("utf-8")
    r2_md5 = calc_md5(args.reads2).decode("utf-8")
    header_qc = "sampleID\ttotReads\tMbases\tavgQual\tspecies\tfractionSpecies\tcheckContamination\tcoverage\tcheckCov\n"
    str_qc = sampleName + "\t" + str(r1["totReads"] + r2["totReads"]) + "\t" + str(
        r1["mbases"] + r2["mbases"]) + "\t" + str(
        round(average(r1["avg_qual"] + r2["avg_qual"]), 2)) + "\t" + sampleSpecies + "\t" + str(round(fraction,3)) + "\t" + checkGenus + "\t" + str(
        int(coverage)) + "X\t" + checkCov + "\n"

    if checkCov == "FAIL" or checkGenus == "FAIL":
        with open("fail_QC.tsv", 'w') as failQC:
            failQC.write(header_qc)
            failQC.write(str_qc)

    str_md5 = "R1_md5\t" + r1_md5 + "\n" + "R2_md5\t" + r2_md5 + "\n"
    print(header_qc + str_qc + "\n" + str_md5)
