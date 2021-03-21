
#!/usr/bin/env python
# coding: utf-8

__author__ = 'LIN Dan'

import argparse
import os
from Bio import SeqIO

# Change format
def turn_bioseq_into_list(file):

    bio_seq = SeqIO.parse(open(file), 'fasta')
    seq_list = list()
    for seq in bio_seq:
        seq_list.append(seq)
    return seq_list

# Count the number of all kinds of gene families
def count_gene_families(seq_list, i, f_single_copy_9_species, f_present_9_species, f_single_copy_7_species, f_present_7_species, f_single_copy_5_species, f_present_5_species, f_single_copy_2_species, f_present_2_species, f_unshared_gene_families, f_singletons_laboriosa, f_singletons_dorsata):

    Temnothorax = list()
    Dinoponera = list()
    Mellifera = list()
    Cerana = list()
    Florea = list()
    Terrestris = list()
    Impatiens = list()
    Laboriosa = list()
    Dorsata = list()

    for seq in seq_list:
        if  'Temnothorax' in seq.id:
            Temnothorax.append(seq)
        if  'Dinoponera' in seq.id:
            Dinoponera.append(seq)
        if  'Mellifera' in seq.id:
            Mellifera.append(seq)
        if  'Cerana' in seq.id:
            Cerana.append(seq)
        if  'Florea' in seq.id:
            Florea.append(seq)
        if  'Terrestris' in seq.id:
            Terrestris.append(seq)
        if  'Impatiens' in seq.id:
            Impatiens.append(seq)
        if  'Laboriosa' in seq.id:
            Laboriosa.append(seq)
        if  'Dorsata' in seq.id:
            Dorsata.append(seq)
        else:
            pass

    if all([len(Temnothorax) == 1, len(Dinoponera) == 1, len(Mellifera) == 1, len(Cerana) == 1, len(Florea) == 1, len(Dorsata) == 1, len(Terrestris) == 1, len(Impatiens) == 1, len(Laboriosa) == 1]):
        f_single_copy_9_species.write(str(i) + '\n')

    if all([len(Temnothorax) > 1, len(Dinoponera) > 1, len(Mellifera) > 1, len(Cerana) > 1, len(Florea) > 1, len(Dorsata) > 1, len(Terrestris) > 1, len(Impatiens) > 1, len(Laboriosa) > 1]):
        f_present_9_species.write(str(i) + '\t' + str(len(Laboriosa)) + '\t' + str(len(Dorsata)) + '\n')

    if all([len(Temnothorax) == 0, len(Dinoponera) == 0, len(Mellifera) == 1, len(Cerana) == 1, len(Florea) == 1, len(Dorsata) == 1, len(Terrestris) == 1, len(Impatiens) == 1, len(Laboriosa) == 1]):
        f_single_copy_7_species.write(str(i))
        f_single_copy_7_species.write('\n')

    if all([len(Temnothorax) == 0, len(Dinoponera) == 0, len(Mellifera) > 1, len(Cerana) > 1, len(Florea) > 1, len(Dorsata) > 1, len(Terrestris) > 1, len(Impatiens) > 1, len(Laboriosa) > 1]):
        f_present_7_species.write(str(i) + '\t' + str(len(Laboriosa)) + '\t' + str(len(Dorsata)) + '\n')

    if all([len(Temnothorax) == 0, len(Dinoponera) == 0, len(Mellifera) == 1, len(Cerana) == 1, len(Florea) == 1, len(Dorsata) == 1, len(Terrestris) == 0, len(Impatiens) == 0, len(Laboriosa) == 1]):
        f_single_copy_5_species.write(str(i) + '\n')

    if all([len(Temnothorax) == 0, len(Dinoponera) == 0, len(Mellifera) > 1, len(Cerana) > 1, len(Florea) > 1, len(Dorsata) > 1, len(Terrestris) == 0, len(Impatiens) == 0, len(Laboriosa) > 1]):
        f_present_5_species.write(str(i) + '\t' + str(len(Laboriosa)) + '\t' + str(len(Dorsata)) + '\n')

    if all([len(Temnothorax) == 0, len(Dinoponera) == 0, len(Mellifera) == 0, len(Cerana) == 0, len(Florea) == 0, len(Dorsata) == 1, len(Terrestris) == 0, len(Impatiens) == 0, len(Laboriosa) == 1]):
        f_single_copy_2_species.write(str(i) + '\n')

    if all([len(Temnothorax) == 0, len(Dinoponera) == 0, len(Mellifera) == 0, len(Cerana) == 0, len(Florea) == 0, len(Dorsata) > 1, len(Terrestris) == 0, len(Impatiens) == 0, len(Laboriosa) > 1]):
        f_present_2_species.write(str(i) + '\t' + str(len(Laboriosa)) + '\t' + str(len(Dorsata)) + '\n')

    if all([len(Temnothorax) + len(Dinoponera) + len(Mellifera) + len(Cerana) + len(Florea) + len(Terrestris) + len(Impatiens) > 0, len(Laboriosa) > 0, len(Dorsata) == 0]):
        f_unshared_gene_families.write(str(i) + '\t' + str(len(Laboriosa)) + '\n')

    if all([len(Temnothorax) + len(Dinoponera) + len(Mellifera) + len(Cerana) + len(Florea) + len(Terrestris) + len(Impatiens) > 0, len(Dorsata) > 0, len(Laboriosa) == 0]):
        f_unshared_gene_families.write(str(i) + '\t' + str(len(Dorsata)) + '\n')

    if all([len(Temnothorax) == 0, len(Dinoponera) == 0, len(Mellifera) == 0, len(Cerana) == 0, len(Florea) == 0, len(Dorsata) == 0, len(Terrestris) == 0, len(Impatiens) == 0, len(Laboriosa) > 0]):
        f_singletons_laboriosa.write(str(i) + '\t' + str(len(Laboriosa)) + '\n')

    if all([len(Temnothorax) == 0, len(Dinoponera) == 0, len(Mellifera) == 0, len(Cerana) == 0, len(Florea) == 0, len(Dorsata) > 0, len(Terrestris) == 0, len(Impatiens) == 0, len(Laboriosa) == 0]):
        f_singletons_dorsata.write(str(i) + '\t' + str(len(Dorsata)) + '\n')

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--gene_family_directory',
                        required = True,
                        type = str,
                        help = 'Path to gene families folder. Please provide the absolute path.')
    parser.add_argument('--file_suffix',
                        required = True,
                        type = str,
                        help = 'Suffix of the gene families files, e.g., cds.muscle')
    parser.add_argument('--gene_family_number',
                        required = True,
                        type = int,
                        help = 'Number of the gene families files')
    args = parser.parse_args()

    gene_family_directory = args.gene_family_directory
    file_suffix = args.file_suffix
    gene_family_number = args.gene_family_number

    os.chdir(gene_family_directory)
    f_single_copy_9_species = open('single_copy_9_species', 'a')
    f_present_9_species = open('present_9_species', 'a')
    f_single_copy_7_species = open('single_copy_7_species', 'a')
    f_present_7_species = open('present_7_species', 'a')
    f_single_copy_5_species = open('single_copy_5_species', 'a')
    f_present_5_species = open('present_5_species', 'a')
    f_single_copy_2_species = open('single_copy_2_species', 'a')
    f_present_2_species = open('present_2_species', 'a')
    f_unshared_gene_families = open('unshared_gene_families', 'a')
    f_singletons_laboriosa = open('singletons_laboriosa', 'a')
    f_singletons_dorsata = open('singletons_dorsata', 'a')

    target_file = [str(i) + file_suffix for i in range(0, gene_family_number)]
    length = len(target_file)
    for i in range(0, length):
        print('Now, openning the %d th gene family' %i)
        seq_list = turn_bioseq_into_list(target_file[i])
        count_gene_families(seq_list, i, f_single_copy_9_species, f_present_9_species, f_single_copy_7_species, f_present_7_species, f_single_copy_5_species, f_present_5_species, f_single_copy_2_species, f_present_2_species, f_unshared_gene_families, f_singletons_laboriosa, f_singletons_dorsata)

    f_single_copy_9_species.close()
    f_single_copy_7_species.close()
    f_single_copy_5_species.close()
    f_single_copy_2_species.close()
    f_present_9_species.close()
    f_present_7_species.close()
    f_present_5_species.close()
    f_present_2_species.close()
    f_unshared_gene_families.close()
    f_singletons_laboriosa.close()
    f_singletons_dorsata.close()

if __name__ == "__main__":
    main()




