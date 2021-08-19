from __future__ import division
import sys
import argparse
import re
import os
import csv
from cyvcf2 import VCF
import subprocess
from collections import OrderedDict
from collections import defaultdict
from collections import defaultdict
import glob
import readline
import pandas as pd
import numpy as np
import timeit
import time
import gc
import datetime
from collections import Counter
from pyfasta import Fasta


# Parse Command line Arguments
parser = argparse.ArgumentParser(
    description='This script will parse Freebayes vcf file and extract allele frequency for LDV positions in LDV_to_Strain_coordinate_map.txt file')
required = parser.add_argument_group('Required arguments')
optional = parser.add_argument_group('Optional arguments')
required.add_argument('-Freebayes_vcf', action='store', dest="Freebayes_vcf",
                      help='Freebayes VCF file')
required.add_argument('-Mutect_vcf', action='store', dest="Mutect_vcf",
                      help='GATK Mutect2 VCF file')
optional.add_argument('-coordinate_map', action='store', dest="coordinate_map",
                      help='Genome coordiante map file generated with Crossmap')
args = parser.parse_args()


def extract_freebayes():

    #print ldv_pos_map
    #print ("LDV Position, Strain Position, Reference Allele, Reference Allele Frequency, ALT Allele, ALT Allele Frequency, REF allele Depth, ALT Allele Frequency")

    LDV_abund_file = "results/freebayes/%s" % os.path.basename(args.Freebayes_vcf.replace('_raw.vcf', '_Freebayes_LDV_abund_frequency.csv'))
    fp = open(LDV_abund_file, 'w+')
    vcf = VCF('%s' % args.Freebayes_vcf)
    fp.write("LDV Position, Strain Position, Reference Allele, Reference Allele Frequency, ALT Allele, ALT Allele Frequency, REF allele Depth, ALT Allele Depth\n")

    freebayes_single_allele_snp_count = 0
    freebayes_multi_allele_snp_count = 0
    freebayes_multi_snp_positions = []
    freebayes_dict = defaultdict(list)
    freebayes_AF_dict = {}

    for v in vcf:
        #print v.POS
        key = v.POS
        if int(v.INFO.get('NUMALT')) == 0 and int(v.INFO.get('DP')) > 0:
            ALT_allele = "."
            Alt_allele_freq = 0
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = 0
            #print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele, Alt_allele_freq, reference_depth, alt_allele_depth))
            #fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele, Alt_allele_freq, reference_depth, alt_allele_depth))

        elif int(v.INFO.get('NUMALT')) == 1 and int(v.INFO.get('DP')) > 0:
            freebayes_single_allele_snp_count += 1
            ALT_allele = v.ALT[0]
            Alt_allele_freq = int(v.INFO.get('AO')) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO'))
            #print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele, Alt_allele_freq, reference_depth, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele, Alt_allele_freq, reference_depth, alt_allele_depth))
            freebayes_dict[v.POS].append(ALT_allele)
            freebayes_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele, Alt_allele_freq, reference_depth, alt_allele_depth)

        elif int(v.INFO.get('NUMALT')) == 2 and int(v.INFO.get('DP')) > 0:
            # print ("NUMALT > 1")
            # print (str(v.INFO.get('AO')[1]))
            freebayes_multi_allele_snp_count += 1
            freebayes_multi_snp_positions.append(v.POS)
            # Calculate Allele Frequency of first ALT allele
            ALT_allele_1 = v.ALT[0]
            Alt_allele_freq_1 = int(v.INFO.get('AO')[0]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[0])
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth, alt_allele_depth))
            freebayes_dict[v.POS].append(ALT_allele_1)
            freebayes_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth,
            alt_allele_depth)

            # Calculate Allele Frequency of second ALT allele
            ALT_allele_2 = v.ALT[1]
            Alt_allele_freq_2 = int(v.INFO.get('AO')[1]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[1])
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth, alt_allele_depth))
            freebayes_dict[v.POS].append(ALT_allele_2)
            freebayes_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth,
            alt_allele_depth)

        elif int(v.INFO.get('NUMALT')) == 3 and int(v.INFO.get('DP')) > 0:
            freebayes_multi_allele_snp_count += 1
            freebayes_multi_snp_positions.append(v.POS)
            # Calculate Allele Frequency of first ALT allele
            ALT_allele_1 = v.ALT[0]
            Alt_allele_freq_1 = int(v.INFO.get('AO')[0]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[0])
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth, alt_allele_depth))
            freebayes_dict[v.POS].append(ALT_allele_1)
            freebayes_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth,
            alt_allele_depth)

            # Calculate Allele Frequency of second ALT allele
            ALT_allele_2 = v.ALT[1]
            Alt_allele_freq_2 = int(v.INFO.get('AO')[1]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[1])
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth, alt_allele_depth))
            freebayes_dict[v.POS].append(ALT_allele_2)
            freebayes_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth,
            alt_allele_depth)

            # Calculate Allele Frequency of third ALT allele
            ALT_allele_3 = v.ALT[2]
            Alt_allele_freq_3 = int(v.INFO.get('AO')[2]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[2])
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth, alt_allele_depth))
            freebayes_dict[v.POS].append(ALT_allele_3)
            freebayes_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth,
            alt_allele_depth)

        elif int(v.INFO.get('NUMALT')) == 4 and int(v.INFO.get('DP')) > 0:
            freebayes_multi_allele_snp_count += 1
            freebayes_multi_snp_positions.append(v.POS)
            # Calculate Allele Frequency of first ALT allele
            ALT_allele_1 = v.ALT[0]
            Alt_allele_freq_1 = int(v.INFO.get('AO')[0]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[0])
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth, alt_allele_depth))
            freebayes_dict[v.POS].append(ALT_allele_1)
            freebayes_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth,
            alt_allele_depth)

            # Calculate Allele Frequency of second ALT allele
            ALT_allele_2 = v.ALT[1]
            Alt_allele_freq_2 = int(v.INFO.get('AO')[1]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[1])
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth, alt_allele_depth))
            freebayes_dict[v.POS].append(ALT_allele_2)
            freebayes_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth,
            alt_allele_depth)

            # Calculate Allele Frequency of third ALT allele
            ALT_allele_3 = v.ALT[2]
            Alt_allele_freq_3 = int(v.INFO.get('AO')[2]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[2])
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth, alt_allele_depth))
            freebayes_dict[v.POS].append(ALT_allele_3)
            freebayes_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth,
            alt_allele_depth)

            # Calculate Allele Frequency of third ALT allele
            ALT_allele_4 = v.ALT[3]
            Alt_allele_freq_4 = int(v.INFO.get('AO')[3]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[3])
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_4, Alt_allele_freq_4, reference_depth, alt_allele_depth))
            freebayes_dict[v.POS].append(ALT_allele_4)
            freebayes_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_4, Alt_allele_freq_4, reference_depth,
            alt_allele_depth)

        elif int(v.INFO.get('NUMALT')) == 5 and int(v.INFO.get('DP')) > 0:
            freebayes_multi_allele_snp_count += 1
            freebayes_multi_snp_positions.append(v.POS)
            # Calculate Allele Frequency of first ALT allele
            ALT_allele_1 = v.ALT[0]
            Alt_allele_freq_1 = int(v.INFO.get('AO')[0]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[0])
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth, alt_allele_depth))
            freebayes_dict[v.POS].append(ALT_allele_1)
            freebayes_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth,
            alt_allele_depth)

            # Calculate Allele Frequency of second ALT allele
            ALT_allele_2 = v.ALT[1]
            Alt_allele_freq_2 = int(v.INFO.get('AO')[1]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[1])
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth, alt_allele_depth))
            freebayes_dict[v.POS].append(ALT_allele_2)
            freebayes_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth,
            alt_allele_depth)

            # Calculate Allele Frequency of third ALT allele
            ALT_allele_3 = v.ALT[2]
            Alt_allele_freq_3 = int(v.INFO.get('AO')[2]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[2])
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth, alt_allele_depth))
            freebayes_dict[v.POS].append(ALT_allele_3)
            freebayes_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth,
            alt_allele_depth)

            # Calculate Allele Frequency of third ALT allele
            ALT_allele_4 = v.ALT[3]
            Alt_allele_freq_4 = int(v.INFO.get('AO')[3]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[3])
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_4, Alt_allele_freq_4, reference_depth, alt_allele_depth))
            freebayes_dict[v.POS].append(ALT_allele_4)
            freebayes_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_4, Alt_allele_freq_4, reference_depth,
            alt_allele_depth)

            # Calculate Allele Frequency of third ALT allele
            ALT_allele_5 = v.ALT[4]
            Alt_allele_freq_5 = int(v.INFO.get('AO')[4]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[4])
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            # key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_5, Alt_allele_freq_5, reference_depth, alt_allele_depth))
            freebayes_dict[v.POS].append(ALT_allele_5)
            freebayes_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_5, Alt_allele_freq_5, reference_depth,
            alt_allele_depth)

    print "Number of Positions with Single SNP allele in Freebayes vcf - %s" % freebayes_single_allele_snp_count
    print "Number of Positions with Multiple SNP allele in Freebayes vcf - %s" % freebayes_multi_allele_snp_count
    #print sorted(freebayes_multi_snp_positions)
    fp.close()
    return freebayes_dict, freebayes_AF_dict, freebayes_single_allele_snp_count, freebayes_multi_allele_snp_count

def compare_variant_calls():
    reference_genome_LDV_abund = "5261-4089-0-RVRE_Aus0004_with_metagenome_LDV_abund_frequency.csv"
    same_genome_LDV_abund = "5261-4089-0-RVRE_same_assembly_with_metagenome_LDV_abund_frequency.csv"
    reference_genome_LDV_abund_map = dict()
    with open("%s" % reference_genome_LDV_abund, 'rU') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader, None)
        for row in csv_reader:
            if row[4] == ".":
                reference_genome_LDV_abund_map[row[0]] = (row[2])
            else:
                reference_genome_LDV_abund_map[row[0]] = (row[4])

    same_genome_LDV_abund_map = dict()
    count_not_found = 0
    with open("%s" % same_genome_LDV_abund, 'rU') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        next(csv_reader, None)
        for row in csv_reader:
            if row[1] != "Not found(Split/Unmap)":
                same_genome_LDV_abund_map[row[0]] = (row[2])
            elif row[1] == "Not found(Split/Unmap)":
                count_not_found += 1
    print ("Length of reference_genome_LDV_abund_map - %s" % len(reference_genome_LDV_abund_map))
    print ("Length of same_genome_LDV_abund - %s" % len(same_genome_LDV_abund_map))

    print ("Not found - %s" % count_not_found)


    for key in reference_genome_LDV_abund_map.keys():
        if key in same_genome_LDV_abund_map.keys():
            if reference_genome_LDV_abund_map[key] == same_genome_LDV_abund_map[key]:
                #continue
                print ("Match: %s, %s, %s" % (key, reference_genome_LDV_abund_map[key], same_genome_LDV_abund_map[key]))
            else:
                #continue
                if reference_genome_LDV_abund_map[key] == ".":
                    print ("Not a Match: %s, %s (No variant called for Aus), %s" % (key, reference_genome_LDV_abund_map[key], same_genome_LDV_abund_map[key]))
                else:
                    print ("Not a Match: %s, %s, %s" % (key, reference_genome_LDV_abund_map[key], same_genome_LDV_abund_map[key]))

        else:
            print ("Key Not found - %s" % key)

class Genotype(object):
    __slots__ = ('alleles', 'phased')

    def __init__(self, li):
        self.alleles = li[:-1]
        self.phased = li[-1]

    def __str__(self):
        sep = "/|"[int(self.phased)]
        return sep.join("0123."[a] for a in self.alleles)
    __repr__ = __str__

def common_member(a, b):
    a_set = set(a)
    b_set = set(b)

    # check length
    if len(a_set.intersection(b_set)) > 0:
        return (a_set.intersection(b_set))
    else:
        return ("no common elements")

def extract_gatk_Mutect():
    #print ("LDV Position, Strain Position, Reference Allele, Reference Allele Frequency, ALT Allele, ALT Allele Frequency, REF allele Depth, ALT Allele Depth")
    LDV_abund_file = "results/gatk_mutect/%s" % os.path.basename(args.Mutect_vcf.replace('_Mutect_raw.vcf', '_GATK_Mutect2_LDV_abund_frequency.csv'))
    fp = open(LDV_abund_file, 'w+')
    vcf = VCF('%s' % args.Mutect_vcf)
    fp.write(
        "LDV Position, Strain Position, Reference Allele, Reference Allele Frequency, ALT Allele, ALT Allele Frequency, REF allele Depth, ALT Allele Depth\n")
    gatk_single_allele_snp_count = 0
    gatk_multi_allele_snp_count = 0
    gatk_multi_snp_positions = []
    gatk_dict = defaultdict(list)
    gatk_AF_dict = {}
    for v in vcf:
        key = v.POS
        alt = v.ALT


        if len(alt) == 1 and alt[0] == "<NON_REF>":
            ALT_allele = "."
            Alt_allele_freq = 0
            Reference_allele = v.REF
            total_depth = int(v.format('DP'))
            if int(v.format('DP')) == 0:
                Reference_allele_freq = 0
                alt_allele_depth = 0
                reference_depth = 0
            elif int(v.format('DP')) > 0:
                Reference_allele_freq = int(v.format('DP')) / int(v.format('DP'))
                reference_depth = int(v.format('DP'))
                alt_allele_depth = 0
            #print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele, Alt_allele_freq, reference_depth, alt_allele_depth))
            #fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele, Alt_allele_freq, reference_depth, alt_allele_depth))
        elif len(alt) == 2:
            gatk_single_allele_snp_count += 1
            ALT_allele = v.ALT[0]
            Reference_allele = v.REF
            AD_array = (((((str(v.format('AD'))).replace('  ', ' ')).replace('[ ', '[')).replace('[', '')).replace(']', '')).split(' ')

            AD_array = list(filter(None, AD_array))
            # print AD_array
            if int(v.format('DP')) == 0:
                total_depth = int(v.format('DP'))
                Reference_allele_freq = 0
                alt_allele_depth = 0
                reference_depth = 0
                alt_allele_depth = 0
            elif int(v.format('DP')) > 0:
                total_depth = int(v.format('DP'))
                Reference_allele_freq = int(AD_array[0]) / total_depth
                alt_allele_depth = int(AD_array[1])
                Alt_allele_freq = int(AD_array[1]) / total_depth
            #print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele, Alt_allele_freq, reference_depth, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele, Alt_allele_freq, reference_depth, alt_allele_depth))
            gatk_dict[v.POS].append(ALT_allele)
            gatk_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele, Alt_allele_freq, reference_depth, alt_allele_depth)

        elif len(alt) == 3:
            gatk_multi_allele_snp_count += 1
            gatk_multi_snp_positions.append(v.POS)
            if int(v.format('DP')) == 0:
                print "Depth at %s was 0. Please check the VCF file" % v.POS
                exit()
            Reference_allele = v.REF
            AD_array = (((((str(v.format('AD'))).replace('  ', ' ')).replace('[ ', '[')).replace('[', '')).replace(']', '')).split(' ')
            AD_array = list(filter(None, AD_array))
            total_depth = int(v.format('DP'))
            Reference_allele_freq = int(AD_array[0]) / total_depth
            alt_allele_depth_1 = int(AD_array[1])
            Alt_allele_freq_1 = int(AD_array[1]) / total_depth
            ALT_allele_1= v.ALT[0]
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            # key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth, alt_allele_depth_1))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth,alt_allele_depth_1))
            gatk_dict[v.POS].append(ALT_allele_1)
            gatk_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth,
                alt_allele_depth_1)

            alt_allele_depth_2 = int(AD_array[2])
            Alt_allele_freq_2 = int(AD_array[2]) / total_depth
            ALT_allele_2 = v.ALT[1]
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth,alt_allele_depth_2))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth,
            alt_allele_depth_2))
            gatk_dict[v.POS].append(ALT_allele_2)
            gatk_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth,
                alt_allele_depth_2)

        elif len(alt) == 4:
            gatk_multi_allele_snp_count += 1
            gatk_multi_snp_positions.append(v.POS)
            if int(v.format('DP')) == 0:
                print "Depth at %s was 0. Please check the VCF file" % v.POS
                exit()
            Reference_allele = v.REF
            AD_array = (((((str(v.format('AD'))).replace('  ', ' ')).replace('[ ', '[')).replace('[', '')).replace(']', '')).split(' ')
            AD_array = list(filter(None, AD_array))
            total_depth = int(v.format('DP'))
            Reference_allele_freq = int(AD_array[0]) / total_depth
            alt_allele_depth_1 = int(AD_array[1])
            Alt_allele_freq_1 = int(AD_array[1]) / total_depth
            ALT_allele_1= v.ALT[0]
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            # key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth, alt_allele_depth_1))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth,alt_allele_depth_1))
            gatk_dict[v.POS].append(ALT_allele_1)
            gatk_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth,
                alt_allele_depth_1)

            alt_allele_depth_2 = int(AD_array[2])
            Alt_allele_freq_2 = int(AD_array[2]) / total_depth
            ALT_allele_2 = v.ALT[1]
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            #     key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth,
            #     alt_allele_depth_2))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth,
            alt_allele_depth_2))
            gatk_dict[v.POS].append(ALT_allele_2)
            gatk_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth,
                alt_allele_depth_2)

            alt_allele_depth_3 = int(AD_array[3])
            Alt_allele_freq_3 = int(AD_array[3]) / total_depth
            ALT_allele_3 = v.ALT[2]
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            #     key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth,
            #     alt_allele_depth_3))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth,
                alt_allele_depth_3))
            gatk_dict[v.POS].append(ALT_allele_3)
            gatk_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth,
                alt_allele_depth_3)

        elif len(alt) == 5:
            gatk_multi_allele_snp_count += 1
            gatk_multi_snp_positions.append(v.POS)
            if int(v.format('DP')) == 0:
                print "Depth at %s was 0. Please check the VCF file" % v.POS
                exit()
            Reference_allele = v.REF
            AD_array = (((((str(v.format('AD'))).replace('  ', ' ')).replace('[ ', '[')).replace('[', '')).replace(']','')).split(' ')
            AD_array = list(filter(None, AD_array))
            total_depth = int(v.format('DP'))
            Reference_allele_freq = int(AD_array[0]) / total_depth
            alt_allele_depth_1 = int(AD_array[1])
            Alt_allele_freq_1 = int(AD_array[1]) / total_depth
            ALT_allele_1 = v.ALT[0]
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            # key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth, alt_allele_depth_1))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth,
            alt_allele_depth_1))
            gatk_dict[v.POS].append(ALT_allele_1)
            gatk_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth,
                alt_allele_depth_1)

            alt_allele_depth_2 = int(AD_array[2])
            Alt_allele_freq_2 = int(AD_array[2]) / total_depth
            ALT_allele_2 = v.ALT[1]
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            #     key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth,
            #     alt_allele_depth_2))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth,
                alt_allele_depth_2))
            gatk_dict[v.POS].append(ALT_allele_2)
            gatk_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth,
                alt_allele_depth_2)

            alt_allele_depth_3 = int(AD_array[3])
            Alt_allele_freq_3 = int(AD_array[3]) / total_depth
            ALT_allele_3 = v.ALT[2]
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            #     key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth,
            #     alt_allele_depth_3))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth,
                alt_allele_depth_3))
            gatk_dict[v.POS].append(ALT_allele_3)
            gatk_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth,
                alt_allele_depth_3)

            alt_allele_depth_4 = int(AD_array[4])
            Alt_allele_freq_4 = int(AD_array[4]) / total_depth
            ALT_allele_4 = v.ALT[3]
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            #     key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth,
            #     alt_allele_depth_3))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_4, Alt_allele_freq_4, reference_depth,
                alt_allele_depth_4))
            gatk_dict[v.POS].append(ALT_allele_4)
            gatk_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_4, Alt_allele_freq_4, reference_depth,
                alt_allele_depth_4)

        elif len(alt) == 6:
            gatk_multi_allele_snp_count += 1
            gatk_multi_snp_positions.append(v.POS)
            if int(v.format('DP')) == 0:
                print "Depth at %s was 0. Please check the VCF file" % v.POS
                exit()
            Reference_allele = v.REF
            AD_array = (((((str(v.format('AD'))).replace('  ', ' ')).replace('[ ', '[')).replace('[', '')).replace(']',
                                                                                                                   '')).split(
                ' ')
            AD_array = list(filter(None, AD_array))
            total_depth = int(v.format('DP'))
            Reference_allele_freq = int(AD_array[0]) / total_depth
            alt_allele_depth_1 = int(AD_array[1])
            Alt_allele_freq_1 = int(AD_array[1]) / total_depth
            ALT_allele_1 = v.ALT[0]
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            # key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth, alt_allele_depth_1))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth,
                alt_allele_depth_1))
            gatk_dict[v.POS].append(ALT_allele_1)
            gatk_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, reference_depth,
                alt_allele_depth_1)

            alt_allele_depth_2 = int(AD_array[2])
            Alt_allele_freq_2 = int(AD_array[2]) / total_depth
            ALT_allele_2 = v.ALT[1]
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            #     key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth,
            #     alt_allele_depth_2))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth,
                alt_allele_depth_2))
            gatk_dict[v.POS].append(ALT_allele_2)
            gatk_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, reference_depth,
                alt_allele_depth_2)

            alt_allele_depth_3 = int(AD_array[3])
            Alt_allele_freq_3 = int(AD_array[3]) / total_depth
            ALT_allele_3 = v.ALT[2]
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            #     key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth,
            #     alt_allele_depth_3))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth,
                alt_allele_depth_3))
            gatk_dict[v.POS].append(ALT_allele_3)
            gatk_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth,
                alt_allele_depth_3)

            alt_allele_depth_4 = int(AD_array[4])
            Alt_allele_freq_4 = int(AD_array[4]) / total_depth
            ALT_allele_4 = v.ALT[3]
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            #     key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth,
            #     alt_allele_depth_3))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_4, Alt_allele_freq_4, reference_depth,
                alt_allele_depth_4))
            gatk_dict[v.POS].append(ALT_allele_4)
            gatk_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_4, Alt_allele_freq_4, reference_depth,
                alt_allele_depth_4)

            alt_allele_depth_5 = int(AD_array[5])
            Alt_allele_freq_5 = int(AD_array[5]) / total_depth
            ALT_allele_5 = v.ALT[4]
            # print ("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
            #     key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, reference_depth,
            #     alt_allele_depth_3))
            fp.write("%s,%s,%s,%s,%s,%s,%s,%s\n" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_5, Alt_allele_freq_5, reference_depth,
                alt_allele_depth_5))
            gatk_dict[v.POS].append(ALT_allele_5)
            gatk_AF_dict[v.POS] = "%s,%s,%s,%s,%s,%s,%s,%s" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_5, Alt_allele_freq_5, reference_depth,
                alt_allele_depth_5)
        elif len(alt) == 7:
            print "More than 6 alleles found. Please add a rule to the source code"
            exit()
    print "Number of Positions with Single SNP allele in GATK Mutect vcf - %s" % gatk_single_allele_snp_count
    print "Number of Positions with Multiple SNP allele in GATK Mutect vcf - %s" % gatk_multi_allele_snp_count

    #print sorted(gatk_multi_snp_positions)
    return gatk_dict, gatk_AF_dict, gatk_single_allele_snp_count, gatk_multi_allele_snp_count


def generate_Final_LDV_abund(gatk_keys, freebayes_keys, uniq_keys):
    positions_called_by_both = 0
    position_called_by_only_freebayes = 0
    position_called_by_only_gatk = 0
    positions_called_by_both_and_same = 0
    positions_called_by_both_but_different = 0

    Final_LDV_abund_file = "results/MergedVCF/%s" % os.path.basename(
        args.Mutect_vcf.replace('_raw.vcf', '_Final_LDV_abund_frequency.csv'))
    fp = open(Final_LDV_abund_file, 'w+')

    header = "LDV Position, Strain Position, Freebayes Reference Allele, Freebayes Reference Allele Frequency, Freebayes ALT Allele, Freebayes ALT Allele Frequency, Freebayes REF allele Depth, Freebayes ALT Allele Depth, LDV Position, Strain Position, GATK Mutect Reference Allele, GATK Mutect Reference Allele Frequency, GATK Mutect ALT Allele, GATK Mutect ALT Allele Frequency, GATK Mutect REF allele Depth, GATK Mutect ALT Allele Depth, Vote, Comment"

    fp.write("%s\n" % header)
    total_unique_positions = len(uniq_keys)

    for position in sorted(uniq_keys):
        AF_string = ""
        if position in gatk_dict.keys() and position in freebayes_dict.keys():
            positions_called_by_both += 1
            if str(gatk_dict[position]) == str(freebayes_dict[position]):
                positions_called_by_both_and_same += 1
                AF_string = "%s,%s,%s,%s\n" % (
                freebayes_AF_dict[position], gatk_AF_dict[position], "2/2", "Same Allele in Both")
                # print "%s, %s, %s" % (position, str(gatk_dict[position]), str(freebayes_dict[position]))
            else:
                AF_string = "%s,%s,%s,%s\n" % (
                freebayes_AF_dict[position], gatk_AF_dict[position], "2/2", "Different Allele in Both")
                positions_called_by_both_but_different += 1
                # print "%s, %s, %s" % (position, str(gatk_dict[position]), str(freebayes_dict[position]))
        elif position in gatk_dict.keys() and position not in freebayes_dict.keys():
            AF_string = "%s,%s,%s,%s\n" % (
            "NA,NA,NA,NA,NA,NA,NA,NA", gatk_AF_dict[position], "1/2", "Called only in GATK Mutect2")
            position_called_by_only_gatk += 1

        elif position not in gatk_dict.keys() and position in freebayes_dict.keys():
            AF_string = "%s,%s,%s,%s\n" % (
            freebayes_AF_dict[position], "NA,NA,NA,NA,NA,NA,NA,NA", "1/2", "Called only in Freebayes")
            position_called_by_only_freebayes += 1
        fp.write("%s\n" % AF_string)

    Final_LDV_abund_stats_file = "results/MergedVCF/%s" % os.path.basename(
        args.Mutect_vcf.replace('_Mutect_raw.vcf', '_Final_LDV_abund_frequency_stats.csv'))
    fp2 = open(Final_LDV_abund_stats_file, 'w+')

    header = "total_unique_positions, positions_called_by_both, positions_called_by_both_and_same, positions_called_by_both_but_different, position_called_by_only_freebayes, position_called_by_only_gatk, Number of Positions with Single SNP allele in GATK Mutect vcf, Number of Positions with Multiple SNP allele in GATK Mutect vcf, Number of Positions with Single SNP allele in Freebayes vcf, Number of Positions with Multiple SNP allele in Freebayes vcf\n"

    fp2.write(header)
    fp2.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" % (
    total_unique_positions, positions_called_by_both, positions_called_by_both_and_same,
    positions_called_by_both_but_different, position_called_by_only_freebayes, position_called_by_only_gatk,
    gatk_single_allele_snp_count, gatk_multi_allele_snp_count, freebayes_single_allele_snp_count,
    freebayes_multi_allele_snp_count))

    print "total_unique_positions - %s" % total_unique_positions
    print "positions_called_by_both - %s" % positions_called_by_both
    print "positions_called_by_both_and_same - %s" % positions_called_by_both_and_same
    print "positions_called_by_both_but_different - %s" % positions_called_by_both_but_different
    print "position_called_by_only_freebayes - %s" % position_called_by_only_freebayes
    print "position_called_by_only_gatk - %s" % position_called_by_only_gatk

    fp2.close()

gatk_dict, gatk_AF_dict, gatk_single_allele_snp_count, gatk_multi_allele_snp_count  = extract_gatk_Mutect()

freebayes_dict, freebayes_AF_dict, freebayes_single_allele_snp_count, freebayes_multi_allele_snp_count = extract_freebayes()

gatk_keys = list(gatk_dict.keys())
freebayes_keys = list(freebayes_dict.keys())
uniq_keys = gatk_keys + freebayes_keys

generate_Final_LDV_abund(gatk_keys, freebayes_keys, uniq_keys)