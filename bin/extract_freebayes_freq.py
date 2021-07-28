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
required.add_argument('-vcf', action='store', dest="vcf",
                      help='Freebayes VCF file')
optional.add_argument('-coordinate_map', action='store', dest="coordinate_map",
                      help='Genome coordiante map file generated with Crossmap')
args = parser.parse_args()


def extract_freebayes():
    # ldv_pos_map = dict()
    # with open("%s" % args.coordinate_map, 'rU') as csv_file:
    #     csv_reader = csv.reader(csv_file, delimiter=',')
    #     next(csv_reader, None)
    #     for row in csv_reader:
    #         ldv_pos_map[row[0]] = (row[1]).lstrip()

    #print ldv_pos_map
    print ("LDV Position, Strain Position, Reference Allele, Reference Allele Frequency, ALT Allele, ALT Allele Frequency")

    LDV_abund_file = "%s" % args.vcf.replace('_raw.vcf', '_LDV_abund_frequency.csv')
    fp = open(LDV_abund_file, 'w+')
    vcf = VCF('%s' % args.vcf)
    fp.write("LDV Position, Strain Position, Reference Allele, Reference Allele Frequency, ALT Allele, ALT Allele Frequency, ALT Allele Depth\n")
    # for key, value in ldv_pos_map.iteritems():
    #     if value == "Not found":
    #         print "%s,%s(Split/Unmap),.,.,.,." % (key,value)
    #         fp.write("%s,%s(Split/Unmap),.,.,.,." % (key,value))
    #     else:
    #         #Aus0004 Genome name - c
    #         #Assembly fasta header name - 5261-4089-0-RVRE_assembly

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
            print ("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele, Alt_allele_freq, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele, Alt_allele_freq, alt_allele_depth))

        elif int(v.INFO.get('NUMALT')) == 1 and int(v.INFO.get('DP')) > 0:
            ALT_allele = v.ALT[0]
            Alt_allele_freq = int(v.INFO.get('AO')) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO'))
            print ("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele, Alt_allele_freq, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele, Alt_allele_freq, alt_allele_depth))

        elif int(v.INFO.get('NUMALT')) == 2 and int(v.INFO.get('DP')) > 0:
            print ("NUMALT > 1")
            print (str(v.INFO.get('AO')[1]))

            # Calculate Allele Frequency of first ALT allele
            ALT_allele_1 = v.ALT[0]
            Alt_allele_freq_1 = int(v.INFO.get('AO')[0]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[0])
            print ("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, alt_allele_depth))

            # Calculate Allele Frequency of second ALT allele
            ALT_allele_2 = v.ALT[1]
            Alt_allele_freq_2 = int(v.INFO.get('AO')[1]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[1])
            print ("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, alt_allele_depth))

        elif int(v.INFO.get('NUMALT')) == 3 and int(v.INFO.get('DP')) > 0:

            # Calculate Allele Frequency of first ALT allele
            ALT_allele_1 = v.ALT[0]
            Alt_allele_freq_1 = int(v.INFO.get('AO')[0]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[0])
            print ("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s\n" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, alt_allele_depth))

            # Calculate Allele Frequency of second ALT allele
            ALT_allele_2 = v.ALT[1]
            Alt_allele_freq_2 = int(v.INFO.get('AO')[1]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[1])
            print ("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, alt_allele_depth))

            # Calculate Allele Frequency of third ALT allele
            ALT_allele_3 = v.ALT[2]
            Alt_allele_freq_3 = int(v.INFO.get('AO')[2]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[2])
            print ("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, alt_allele_depth))

        elif int(v.INFO.get('NUMALT')) == 4 and int(v.INFO.get('DP')) > 0:

            # Calculate Allele Frequency of first ALT allele
            ALT_allele_1 = v.ALT[0]
            Alt_allele_freq_1 = int(v.INFO.get('AO')[0]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[0])
            print ("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s\n" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, alt_allele_depth))

            # Calculate Allele Frequency of second ALT allele
            ALT_allele_2 = v.ALT[1]
            Alt_allele_freq_2 = int(v.INFO.get('AO')[1]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[1])
            print ("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, alt_allele_depth))

            # Calculate Allele Frequency of third ALT allele
            ALT_allele_3 = v.ALT[2]
            Alt_allele_freq_3 = int(v.INFO.get('AO')[2]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[2])
            print ("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, alt_allele_depth))

            # Calculate Allele Frequency of third ALT allele
            ALT_allele_4 = v.ALT[3]
            Alt_allele_freq_4 = int(v.INFO.get('AO')[3]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[3])
            print ("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s\n" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_4, Alt_allele_freq_4, alt_allele_depth))

        elif int(v.INFO.get('NUMALT')) == 5 and int(v.INFO.get('DP')) > 0:

            # Calculate Allele Frequency of first ALT allele
            ALT_allele_1 = v.ALT[0]
            Alt_allele_freq_1 = int(v.INFO.get('AO')[0]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[0])
            print ("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s\n" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_1, Alt_allele_freq_1, alt_allele_depth))

            # Calculate Allele Frequency of second ALT allele
            ALT_allele_2 = v.ALT[1]
            Alt_allele_freq_2 = int(v.INFO.get('AO')[1]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[1])
            print ("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_2, Alt_allele_freq_2, alt_allele_depth))

            # Calculate Allele Frequency of third ALT allele
            ALT_allele_3 = v.ALT[2]
            Alt_allele_freq_3 = int(v.INFO.get('AO')[2]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[2])
            print ("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, alt_allele_depth))

            # Calculate Allele Frequency of third ALT allele
            ALT_allele_4 = v.ALT[3]
            Alt_allele_freq_4 = int(v.INFO.get('AO')[3]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[3])
            print ("%s,%s,%s,%s,%s,%s,%s\n" % (key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s\n" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_4, Alt_allele_freq_4, alt_allele_depth))

            # Calculate Allele Frequency of third ALT allele
            ALT_allele_5 = v.ALT[4]
            Alt_allele_freq_5 = int(v.INFO.get('AO')[4]) / int(v.INFO.get('DP'))
            Reference_allele = v.REF
            Reference_allele_freq = int(v.INFO.get('RO')) / int(v.INFO.get('DP'))
            total_depth = int(v.INFO.get('DP'))
            reference_depth = int(v.INFO.get('RO'))
            alt_allele_depth = int(v.INFO.get('AO')[4])
            print ("%s,%s,%s,%s,%s,%s,%s\n" % (
            key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_3, Alt_allele_freq_3, alt_allele_depth))
            fp.write("%s,%s,%s,%s,%s,%s,%s\n" % (
                key, v.POS, Reference_allele, Reference_allele_freq, ALT_allele_5, Alt_allele_freq_5, alt_allele_depth))

    fp.close()

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

extract_freebayes()

#compare_variant_calls()
