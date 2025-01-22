#!/usr/bin/env python3
from __future__ import division
import argparse
import subprocess
import re

# code updated and adapted from https://github.com/kseniaarkhipova/RedRed/blob/master/RedRed

# Below the path to mummer executables can be specified if it is not in PATH variable:
path_to_mummer = ''

##############################################################

parser = argparse.ArgumentParser(description='Script for reducing redundancy and clustering of DNA sequences')
parser.add_argument('--fna', dest='seqs', required=True, nargs='?', help='Multifasta file with DNA sequences (required)')
parser.add_argument('-b', dest='batch', default=-1, nargs='?', type=int, help='Number of sequences to be aligned in each iteration. Setting this parameter too low will reduce efficiency of the script, too high - crash of script due to not enough memory issue. You can start from 10000. Default = all')
parser.add_argument('--pident', dest='identity', default=95, nargs='?', type=float, help='Percentage of nucleotide identity of alignment to consider two sequences as the same (in form 0-100). Default = 95')
parser.add_argument('--coverage', dest='cover', default=80, nargs='?', type=float, help='Part/percentge of shortest of two aligned sequences which should be covered by the alignment (in form 0-100). Default = 80')
parser.add_argument('--lenlim', dest='length_limit', default=1, nargs='?', help='Filter out sequences shorter than the limit. Default = no filtering')
args = parser.parse_args()

prefix = re.search('(\S+)\.', args.seqs).group(1)

sequences = {}
def processing(name, seqs):
    seq = ''.join(seqs)
    if len(seq) >= int(args.length_limit):
        sequences[name] = seq

print('\n\nScript was launched with the following parameters: \n Fasta file name: ' + args.seqs + '\n Number of sequences in each iteration: ' + str(args.batch) + '\n Percentage of nucl. identity: ' + str(args.identity) + '\n Coverage of shortest sequence: ' + str(args.cover) + '\n Sequences shorter than ' + str(args.length_limit) + ' bp will be filtered out \n\nPresorting is starting...\n\n')

fl = 0
with open(args.seqs) as file1:
    for line in file1:
        if line.startswith('>') and fl == 0:
            l = []
            name = re.search('>(.*)$', line.strip()).group(1)
            fl = 1
        elif line.startswith('>') and fl == 1:
            processing(name, l)
            l = []
            name = re.search('>(\S+)', line.strip()).group(1)
        elif not line.startswith('>'):
            l.append(line.strip())
    processing(name, l)

p = sorted(sequences, key=lambda record:len(sequences[record]), reverse=True)

with open('temporal_1.fasta', 'w') as file2:
    for item in p:
        file2.write('>' + item + '\n' + sequences[item] + '\n')

del sequences
seen_dict = {}
non_red = []
seen = set()

f = open(prefix + '_clusters.txt', 'w')
f.close()

all_contigs = subprocess.check_output(['grep', '>', 'temporal_1.fasta'], universal_newlines=True)
num_cont = len(all_contigs.split())

print('Presorting has finished. Iterative alignment is starting...\n\n')
if args.batch == -1:
    args.batch = num_cont
count = 1
while len(seen_dict) != num_cont:
    print('Iteration: ' + str(count) + '; Number of analysed sequences: ' + str(len(seen_dict)) + '; Total number of sequences: ' + str(num_cont) + '\n\n')
    seqs_list = []
    name = ''
    with open("temporal_1.fasta") as file1:
        for line in file1:
            if line.startswith('>') and re.search('>(.*)$', line).group(1) not in seen_dict:
                if l and name:
                    d = {}
                    d[name] = ''.join(l)
                    seqs_list.append(d)
                l = []
                name = re.search('>(.*)$', line.strip()).group(1)
                flag = 1
            elif line.startswith('>') and re.search('>(.*)$', line).group(1) in seen_dict:
                if l and name:
                    d = {}
                    d[name] = ''.join(l)
                    seqs_list.append(d)
                l = []
                flag = 0
            elif not line.startswith('>') and flag == 1:
                l.append(line.strip())
    if l and name:
        d = {}
        d[name] = ''.join(l)
        seqs_list.append(d)
    count2 = 1
    query_out = open('temporal_1.fasta', 'w')
    ref_out = open('temporal_2.fasta', 'w')
    left = []
    for item in seqs_list:
        # Python 3 update: dict.keys() returns a view object, convert to list
        left.append(list(item.keys())[0])
        query_out.write('>' + list(item.keys())[0] + '\n' + list(item.values())[0] + '\n')
        if count2 <= int(args.batch):
            ref_out.write('>' + list(item.keys())[0] + '\n' + list(item.values())[0] + list(item.values())[0] + '\n')
            count2 += 1
    query_out.close()
    ref_out.close()
    del seqs_list

    print('Alignment with mummer\n')
    comm1 = (path_to_mummer + 'nucmer -maxmatch --nooptimize -p nucmer temporal_2.fasta temporal_1.fasta').split()
    comm2 = (path_to_mummer + 'show-coords -r -l -T -H nucmer.delta').split()

    subprocess.call(comm1)
    subprocess.call(comm2, stdout=open('nucmer.coords','w'))    

    print('\nParsing of the alignment result...\n')
    subprocess.call('sort -k8,8rn nucmer.coords > nucmer_sorted.coords', shell=True)
    with open('nucmer_sorted.coords') as parse_file:
        main_dict = {}
        check = set()
        for line in parse_file:
            dat = line.strip().split('\t')
            if dat[9] != dat[10] and float(dat[6]) > int(args.identity) and (int(dat[5]) / int(dat[8]) * 100) > int(args.cover) and (dat[9] not in main_dict) and (dat[9] not in check):
                main_dict[dat[9]] = [dat[10]]
                check.add(dat[10])
                check.add(dat[9])
            elif dat[9] != dat[10] and float(dat[6]) > int(args.identity) and (int(dat[5]) / int(dat[8]) * 100) > int(args.cover) and dat[9] in main_dict:
                main_dict[dat[9]].append(dat[10])
                check.add(dat[10])
    names_ref = left[:int(args.batch)]
    total = []

    while names_ref:
        if names_ref[0] in main_dict:
            clust_set = list(set(main_dict[names_ref[0]]))    
            clust_set.insert(0, names_ref[0])    
            total.append(clust_set)
            non_red.append(names_ref[0])
            seen.update(clust_set)
        elif names_ref[0] not in main_dict and names_ref[0] in check:
            seen.add(names_ref[0])
        else:
            non_red.append(names_ref[0])
            total.append([names_ref[0]])
            seen.add(names_ref[0])
        seen_dict = {z: 1 for z in seen}
        names_ref = [i for i in names_ref if i not in seen_dict]
    with open(prefix + '_clusters.txt','a') as file1:
        for item in total:
            file1.write('\t'.join(item) + '\n')
    count += 1

print('Almost done!')

name_duct = {('>' + z): 1 for z in non_red}
total = []
with open(args.seqs) as file1:
    for line in file1:
        if line.startswith('>') and line.strip() in name_duct:
            total.append(line)
            flag = 1
        elif line.startswith('>') and line.strip() not in name_duct:
            flag = 0
        elif not line.startswith('>') and flag == 1:
            total.append(line)
with open(prefix + "_non_redundant.fasta", 'w') as file2:
    for item in total:
        file2.write(item)

subprocess.call(['rm', '-f', 'nucmer.coords'])
subprocess.call(['rm', '-f', 'nucmer_sorted.coords'])
subprocess.call(['rm', '-f', 'nucmer.delta'])
subprocess.call(['rm', '-f', 'temporal_1.fasta'])
subprocess.call(['rm', '-f', 'temporal_2.fasta'])


print(f"\n\n >>>>> output saved in {prefix}\n\n")