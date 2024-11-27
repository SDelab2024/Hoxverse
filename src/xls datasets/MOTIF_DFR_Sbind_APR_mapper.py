# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 02:04:20 2024

@author: jishn
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 23:51:37 2024

@author: jishn
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Jun  8 17:37:12 2024

@author: jishn
"""

import sys
import pandas as pd
import glob
import os
import csv
import matplotlib.pyplot as plt
import re
from matplotlib.lines import Line2D
import itertools

nes_motifs = [
    r'L.{2,3}[LIVFM].{2,3}L.{1}[LI]',
    r'L.{2,3}[LIVFM].{2,3}[LI]',
    r'L.{2,3}[LIVFM].{2,3}L.{2,3}[LI]',
    r'L.{2,3}[LIVFM].{2,3}L.{2,3}[LIVFM].{2,3}L.{1}[LI]',
    # Add more NES motifs as needed
]      
lig_motifs = [
r'G.{3}[GA]',
r'G.{1}G.{2}[GA]',
r'G.{2}G.{1}[GA]',
r'G.{3}G.{1}',
r'G.{3}GA',
r'G.{4}[GA]',
]
NLS_patterns = {
r'PKKKRKV',
r'KRPAATKKAGQAKKKK',
r'KRK.{0,2}KXK.{0,2}KK',
r'K.{0,2}R.{10,12}K.{0,2}R.{10,12}K.{0,2}K',
r'K.{0,2}R.{0,2}K.{0,2}K',
r'K.{0,2}R.{0,2}P.{0,2}R',
r'K.{0,2}R.{0,2}X{2}K{2}'
}
proline_motifs = {
 r'P.{2}P',
 r'P.{3}P',
 r'PPP',
 r'PP.{2}PP',
 r'P.{2}P.{1}'
}
sh2_patterns = [
r'[YFW][YFWVLI][LI].{2,3}[DE]', # Core SH2 domain consensus sequence
r'pY.{2}P' # Phosphorylated tyrosine motif
]
sh3_patterns = [
r'P.{2}P'
    ]
zinc_finger_motif = [
r'C.{2,4}C.{12}H.{3,5}H',   
    ]
leucine_zipper_motif = [
r'L.{6}L.{6}L.{6}L',
    ]
Walker_A_motif = [
r'G.{4}GK[TS]'   
    ]
Walker_B_motif = [
r'[AFILMVWY]{4}DE'   
    ]
pdz_pattern = r'[A-Z]-[ST]-[A-ZVL]'

common_MOTIFs = [nes_motifs, lig_motifs, NLS_patterns, pdz_pattern, proline_motifs, sh2_patterns, sh3_patterns, zinc_finger_motif, leucine_zipper_motif, Walker_A_motif, Walker_B_motif, pdz_pattern]
MOTIF_headings = ['NES Motif', 'LIG Motif', 'NLS Motif', 'PDZ Motif', 'Proline Motif', 'SH2 Motif', 'SH3 Motif', 'Zinc Finger Motif', 'Leucine Zipper Motif', 'Walker A Motif', 'Walker B Motif']


def APR_mapper(APR_filename, value_cuttoff, value_repeat_num):
        def create_APR_map():
            return {
                'protein_names': [],
                'residue_lists': [],
                'amino_acid_lists': [],
                'Values_lists': [],
                'weight_lists': [],
                'sequence_length': []
            }
        def add_APR_map(data, protein_name, residues, amino_acids, APR_values, weights, seq_len):
            data['protein_names'].append(protein_name)
            data['residue_lists'].append(residues)
            data['amino_acid_lists'].append(amino_acids)
            data['Values_lists'].append(APR_values)
            data['weight_lists'].append(weights)
            data['sequence_length'].append(seq_len)

        APR_map = create_APR_map()
        with open(APR_filename, 'r') as file:
            lines = file.readlines()
        num_seqs = 0
        for line in lines:   
            num_seqs += line.count(">")
        temp = num_seqs
        for line in lines:   
            if line.startswith('>'):
                # Extract the protein name
                protein_name = line.split('|')[2]
                end = protein_name.find("_")
                protein_name = protein_name[0:end]
                residue_count = 0
                APR_residue_list = []
                APR_weight_list = []
                APR_values = []
                APR_amino_acid_list =[]
                key = 0
                for next_line in lines[lines.index(line) + 1:]:
                    if next_line.startswith('>'):
                        break
                    if next_line.count("Num") > 0 or next_line.count("-") > 0:
                        continue
                    else:
                        residue_count += 1
                        if next_line.count("f"):
                            temp_next_line = ""
                            for k in next_line:
                                if k != "f":
                                    temp_next_line += k
                                else:
                                    temp_next_line += " "
                            next_line = temp_next_line
                        columns = next_line.split()
                        value = float(columns[2])
                        amino_acid = columns[1]
                        if value >= value_cuttoff:
                            key += 1
                            APR_amino_acid_list.append(amino_acid)
                            APR_residue_list.append(residue_count)
                            APR_weight_list.append(temp)
                            APR_values.append(value)
                        else:
                            if key < value_repeat_num and key > 0:
                                list_len = len(APR_amino_acid_list)
                                APR_amino_acid_list = APR_amino_acid_list[0:list_len-key]
                                APR_residue_list = APR_residue_list[0:list_len-key]
                                APR_weight_list = APR_weight_list[0:list_len-key]
                                APR_values = APR_values[0:list_len-key]
                            key = 0
                seq_len = residue_count
                add_APR_map(APR_map, protein_name, APR_residue_list, APR_amino_acid_list, APR_values, APR_weight_list, seq_len)
                temp -= 1
        return(APR_map)
                            
def Sbind_DFR_mapper(seq_filename, Sbind_DFR_filename, sbind_value_cuttoff, sbind_value_repeat_num, pDP_value_cuttoff):
        with open(Sbind_DFR_filename, 'r') as file:
            lines = file.readlines()
        num_seqs = 0
        def create_Sbind_DFR_map():
            return {
                'protein_names': [],
                'residue_lists': [],
                'amino_acid_lists': [],
                'pDP_values_lists': [],
                'Sbind_values_lists': [],
                'weight_lists': [],
                'sequence_length': []
            }
        def add_Sbind_DFR_map(data, protein_name, residues, amino_acids, pDP_values, Sbind_values, weights, seq_len):
            data['protein_names'].append(protein_name)
            data['residue_lists'].append(residues)
            data['amino_acid_lists'].append(amino_acids)
            data['pDP_values_lists'].append(pDP_values)
            data['Sbind_values_lists'].append(Sbind_values)
            data['weight_lists'].append(weights)
            data['sequence_length'].append(seq_len)
        Sbind_map = create_Sbind_DFR_map()
        DFR_map = create_Sbind_DFR_map()
        for line in lines:   
            num_seqs += line.count(">")
        temp = num_seqs
        for line in lines:   
            if line.startswith('>'):
                # Extract the protein name
                protein_name = line.split('|')[2]
                end = protein_name.find("_")
                protein_name = protein_name[0:end]
                residue_count = 0
                DFR_residue_list = []
                DFR_weight_list = []
                DFR_pDP_values = []
                DFR_Sbind_values = []
                DFR_amino_acid_list =[]
                Sbind_residue_list = []
                Sbind_weight_list = []
                Sbind_pDP_values = []
                Sbind_Sbind_values = []
                Sbind_amino_acid_list =[]
                key = 0
                for next_line in lines[lines.index(line) + 1:]:
                    if next_line.startswith('>'):
                        break
                    if next_line.startswith("Loc") or next_line.startswith(">"):
                        continue
                    else:
                        residue_count += 1
                        columns = next_line.split()
                        # Extract the pDP value (3rd column) and convert to float
                        pDP_value = float(columns[2])
                        sbind_value = float(columns[3])
                        amino_acid = columns[1]
                        if pDP_value >= pDP_value_cuttoff:
                            DFR_amino_acid_list.append(amino_acid)
                            DFR_residue_list.append(residue_count)
                            DFR_weight_list.append(temp)
                            DFR_pDP_values.append(pDP_value)
                        if sbind_value >= sbind_value_cuttoff:
                            key += 1
                            Sbind_amino_acid_list.append(amino_acid)
                            Sbind_residue_list.append(residue_count)
                            Sbind_weight_list.append(temp)
                            Sbind_Sbind_values.append(sbind_value)
                        else:
                            if key < sbind_value_repeat_num and key > 0:
                                list_len = len(Sbind_amino_acid_list)
                                Sbind_amino_acid_list = Sbind_amino_acid_list[0:list_len-key]
                                Sbind_residue_list = Sbind_residue_list[0:list_len-key]
                                Sbind_weight_list = Sbind_weight_list[0:list_len-key]
                                Sbind_Sbind_values = Sbind_Sbind_values[0:list_len-key]
                            key = 0
                seq_len = residue_count                           
                add_Sbind_DFR_map(DFR_map, protein_name, DFR_residue_list, DFR_amino_acid_list, DFR_pDP_values, DFR_Sbind_values, DFR_weight_list, seq_len)
                add_Sbind_DFR_map(Sbind_map, protein_name, Sbind_residue_list, Sbind_amino_acid_list, Sbind_pDP_values, Sbind_Sbind_values, Sbind_weight_list, seq_len)
                temp -= 1
        Sbind_DFR_map = [Sbind_map, DFR_map]
        return(Sbind_DFR_map)
       
def motif_mapper(seq_filename):
    def create_residue_map():
        return {
            'protein_names': [],
            'residue_lists': [],
            'amino_acid_lists': [],
            'weight_lists': [],
            'MOTIF type': [],
            'sequence_length': []

        }
    def add_residue_map(data, protein_name, residues, amino_acids, weights, MOTIF_type, seq_len):
        data['protein_names'].append(protein_name)
        data['residue_lists'].append(residues)
        data['amino_acid_lists'].append(amino_acids)
        data['weight_lists'].append(weights)
        data['MOTIF type'].append(MOTIF_type)
        data['sequence_length'].append(seq_len)

    motif_input = int(input("What motif do you want to search for:\n1.NES motif\n2.LIG motif\n3.NLS motif\n4.PDZ motif\n5.Proline motifs\n6.SH2 Motifs\n7.SH3 Motif\n8.Zinc Finger Motif\n9.Leucine Zipper Motif\n10.Walker A Motif\n11.Walker B Motif\n12.Custom\n"))
    motif_map = create_residue_map()  # Assuming you have a function create_residue_map() defined
    match motif_input:
        case 1:
            motif_of_interest = nes_motifs
            MOTIF_type = 'NES Motif'
        case 2:
            motif_of_interest = lig_motifs
            MOTIF_type = 'LIG Motif'
        case 3:
            motif_of_interest = NLS_patterns
            MOTIF_type = 'NLS Motif'
        case 4:
            motif_of_interest = pdz_pattern
            MOTIF_type = 'PDZ Motif'
        case 5:
            motif_of_interest = proline_motifs
            MOTIF_type = 'Proline Motif'
        case 6:
            motif_of_interest = sh2_patterns
            MOTIF_type = 'SH2 Motif'
        case 7:
            motif_of_interest = sh3_patterns
            MOTIF_type = 'SH3 Motif'
        case 8:
            motif_of_interest = zinc_finger_motif
            MOTIF_type = 'Zinc Finger Motif'
        case 9:
            motif_of_interest = leucine_zipper_motif
            MOTIF_type = 'Leucine Zipper Motif'
        case 10:
            motif_of_interest = Walker_A_motif
            MOTIF_type = 'Walker A Motif'
        case 11:
            motif_of_interest = Walker_B_motif
            MOTIF_type = 'Walker B Motif'
        case 12:
            motif_user_input = input("Enter custom motif: ")
            motif_of_interest = {rf'{motif_user_input}'}
            MOTIF_type = motif_user_input 
    motif_map['MOTIF type'].append(MOTIF_type)
    with open(seq_filename, 'r') as file:
        lines = file.readlines()

    num_seqs = sum(1 for line in lines if line.startswith('>'))

    temp = num_seqs

    for line in lines:
        if line.startswith('>'):
            no_matches_flag = True
            index = 0
            protein_name = line.split('|')[2].split('_')[0]
            residue_list = []
            weight_list = []
            amino_acid_list = []
            for next_line in lines[lines.index(line) + 1:]:
                if next_line.startswith('>'):
                    break
                for char in next_line.strip():
                    index += 1
                    amino_acid_list.append(char)
                    residue_list.append(index)  # Assuming index represents residue number
                    weight_list.append(temp)   # Assuming temp represents some weight
            sequence = ''.join(amino_acid_list)
            seq_len = len(sequence)
            for motif in motif_of_interest:
                matches = re.finditer(motif, sequence)
                for match in matches:
                    no_matches_flag = False
                    start_index = match.start() 
                    end_index = match.end()
                    motif_map['protein_names'].append(protein_name)
                    motif_map['residue_lists'].append(residue_list[start_index:end_index])
                    motif_map['amino_acid_lists'].append(amino_acid_list[start_index:end_index])
                    motif_map['weight_lists'].append(weight_list[start_index:end_index])
                    motif_map['sequence_length'].append(seq_len)
            if no_matches_flag == True:
                motif_map['protein_names'].append(protein_name)
                motif_map['residue_lists'].append([0])
                motif_map['amino_acid_lists'].append([" "])
                motif_map['weight_lists'].append([temp])
                motif_map['sequence_length'].append(seq_len)
            temp -= 1
    
    names = motif_map['protein_names']
    all_res_lists = motif_map['residue_lists']
    all_aa_lists = motif_map['amino_acid_lists']
    all_weight_lists = motif_map['weight_lists']
    all_seq_list = motif_map['sequence_length']
    index = 0
    prev_name = names[index]
    name = prev_name
    motif_map = create_residue_map()
    motif_map['MOTIF type'].append(MOTIF_type)
    while index < len(names):
        res_list = []
        aa_list = []
        weight_list = []
        seq = []
        while name == prev_name:
            for a, b, c in zip(all_res_lists[index], all_aa_lists[index], all_weight_lists[index]):
                res_list.append(a)
                aa_list.append(b)
                weight_list.append(c)
            if index < len(names) -1: 
                curr_name = name
                index += 1
                name = names[index]
            else: 
                name = ''
                index += 1
        seq.append(all_seq_list[index-1])
        unique_residues = []
        duplicate_indices = []
        for pindex, res in enumerate(res_list):
            if res not in unique_residues:
                unique_residues.append(res)
            else:
                duplicate_indices.append(pindex)
        for pindex in sorted(duplicate_indices, reverse=True):
            del res_list[pindex]
            del aa_list[pindex]
            del weight_list[pindex]
        motif_map['protein_names'].append(curr_name)
        motif_map['residue_lists'].append(res_list)
        motif_map['amino_acid_lists'].append(aa_list)
        motif_map['weight_lists'].append(weight_list)
        motif_map['sequence_length'].append(seq)
        prev_name = name
         
    return motif_map
  
def make_CSV(residue_map, Sbind_DFR_map, APR_map, directory, csv_plot_output, DFR_action, Sbind_action, APR_action, MOTIF_action):     
        seq_flag = False
        DFR_in_total_all = []
        Sbind_in_total_all = []
        APR_in_total_all = []
        MOTIF_in_total_all = []
        if MOTIF_action == True: 
            protein_names = residue_map['protein_names']
            residue_lists = residue_map['residue_lists']
            amino_acid_lists = residue_map['amino_acid_lists']
            weight_lists = residue_map['weight_lists']
            data = {'Protein Name': [], 'Residues': [], 'Amino_acids': [], 'Weights': []}
            for protein_name, residues, amino_acids, weights in zip(protein_names, residue_lists, amino_acid_lists, weight_lists):
                data['Protein Name'].extend([protein_name] * len(residues))
                data['Residues'].extend(residues)
                data['Amino_acids'].extend(amino_acids)
                data['Weights'].extend(weights)
            df = pd.DataFrame(data)
            output_file_path = 'new_protein_data.csv'
            df.to_csv(output_file_path, index=False)
            MOTIF_name = residue_map['MOTIF type']
            temp_name = ''
            for i in MOTIF_name[0]:
                temp_name += i
            if MOTIF_headings.count(temp_name) != 0:
                j = MOTIF_headings.index(temp_name)
                MOTIF_defintion = common_MOTIFs[j]
            else:
                MOTIF_defintion = MOTIF_name
                MOTIF_name = temp_name + " Repeats"
            MOTIF_df = pd.DataFrame(MOTIF_defintion, columns=[MOTIF_name])
            output_file_path = 'new_protein_data_MOTIF.csv'
            MOTIF_df.to_csv(output_file_path, index=False)
            sequence_length = residue_map['sequence_length']
            temp_seq_len = []
            for a in sequence_length:
                for b in a:
                    temp_seq_len.append(b)
            sequence_length = temp_seq_len
            print(sequence_length)
            seq_len_df = pd.DataFrame({'Protein Name': protein_names, 'Sequence_Length': sequence_length})        
            seq_flag = True
            count = 0
            while count < len(sequence_length):
                array1 = list(range(1, sequence_length[count] + 1))
                array2 = residue_lists[count]
                common_elements = list(filter(lambda x: x in array2, array1))
                MOTIF_in_total_all.append(len(common_elements)*100/sequence_length[count])
                count += 1 
            protein_names_all_2params = protein_names

        if DFR_action == True: 
            DFR_map = Sbind_DFR_map[1]
            protein_names_DFR = DFR_map['protein_names']
            DFR_residue_lists = DFR_map['residue_lists']
            DFR_amino_acid_lists = DFR_map['amino_acid_lists']
            DFR_weight_lists = DFR_map['weight_lists']
            pDP_values_lists = DFR_map['pDP_values_lists']
            data_DFR = {'Protein Name': [], 'DFR_Residues': [] , 'DFR_Amino_Acids': [], 'pDP_values': [], 'DFR_Weights': []}
            for protein_name, DFR_residues, DFR_amino_acids, pDP_values, DFR_weights in zip(protein_names_DFR, DFR_residue_lists, DFR_amino_acid_lists, pDP_values_lists, DFR_weight_lists):
                data_DFR['Protein Name'].extend([protein_name] * len(DFR_residues))
                data_DFR['DFR_Residues'].extend(DFR_residues)
                data_DFR['DFR_Amino_Acids'].extend(DFR_amino_acids)
                data_DFR['pDP_values'].extend(pDP_values)
                data_DFR['DFR_Weights'].extend(DFR_weights)
            df_DFR = pd.DataFrame(data_DFR)
            output_file_path = 'new_protein_data_DFR.csv'
            df_DFR.to_csv(output_file_path, index=False)
            sequence_length = DFR_map['sequence_length']
            if seq_flag == False:
                seq_len_df = pd.DataFrame({'Protein Name': protein_names_DFR, 'Sequence_Length': sequence_length})        
                seq_flag = True
            count = 0
            while count < len(sequence_length):
                array1 = list(range(1, sequence_length[count] + 1))
                array2 = DFR_residue_lists[count]
                common_elements = list(filter(lambda x: x in array2, array1))
                DFR_in_total_all.append(len(common_elements)*100/sequence_length[count])
                count += 1 
            protein_names_all_2params = protein_names_DFR

        if Sbind_action == True:        
            Sbind_map = Sbind_DFR_map[0]
            protein_names_Sbind = Sbind_map['protein_names']
            Sbind_residue_lists = Sbind_map['residue_lists']
            Sbind_amino_acid_lists = Sbind_map['amino_acid_lists']
            Sbind_weight_lists = Sbind_map['weight_lists']
            Sbind_values_lists = Sbind_map['Sbind_values_lists']
            data_Sbind = {'Protein Name': [], 'Sbind_Residues': [] , 'Sbind_Amino_Acids': [], 'Sbind_values': [], 'Sbind_Weights': []}
            for protein_name, Sbind_residues, Sbind_amino_acids, Sbind_values, Sbind_weights in zip(protein_names_Sbind, Sbind_residue_lists, Sbind_amino_acid_lists, Sbind_values_lists, Sbind_weight_lists):
                data_Sbind['Protein Name'].extend([protein_name] * len(Sbind_residues))
                data_Sbind['Sbind_Residues'].extend(Sbind_residues)
                data_Sbind['Sbind_Amino_Acids'].extend(Sbind_amino_acids)
                data_Sbind['Sbind_values'].extend(Sbind_values)
                data_Sbind['Sbind_Weights'].extend(Sbind_weights)
            df_Sbind = pd.DataFrame(data_Sbind)
            output_file_path = 'new_protein_data_Sbind.csv'
            df_Sbind.to_csv(output_file_path, index=False)
            sequence_length = Sbind_map['sequence_length']
            if seq_flag == False:
                seq_len_df = pd.DataFrame({'Protein Name': protein_names_Sbind, 'Sequence_Length': sequence_length})        
                seq_flag = True
            count = 0
            while count < len(sequence_length):
                array1 = list(range(1, sequence_length[count] + 1))
                array2 = Sbind_residue_lists[count]
                common_elements = list(filter(lambda x: x in array2, array1))
                Sbind_in_total_all.append(len(common_elements)*100/sequence_length[count])
                count += 1   
            protein_names_all_2params = protein_names_Sbind

        if APR_action == True:
            protein_names_APR = APR_map['protein_names']
            APR_residue_lists = APR_map['residue_lists']
            APR_amino_acid_lists = APR_map['amino_acid_lists']
            APR_weight_lists = APR_map['weight_lists']
            APR_values_lists = APR_map['Values_lists']
            data_APR = {'Protein Name': [], 'APR_Residues': [] , 'APR_Amino_Acids': [], 'APR_values': [], 'APR_Weights': []}
            for protein_name, APR_residues, APR_amino_acids, APR_values, APR_weights in zip(protein_names_APR, APR_residue_lists, APR_amino_acid_lists, APR_values_lists, APR_weight_lists):
                data_APR['Protein Name'].extend([protein_name] * len(APR_residues))
                data_APR['APR_Residues'].extend(APR_residues)
                data_APR['APR_Amino_Acids'].extend(APR_amino_acids)
                data_APR['APR_values'].extend(APR_values)
                data_APR['APR_Weights'].extend(APR_weights)
            df_APR = pd.DataFrame(data_APR)
            output_file_path = 'new_protein_data_APR.csv'
            df_APR.to_csv(output_file_path, index=False)
            sequence_length = APR_map['sequence_length']
            if seq_flag == False:
                seq_len_df = pd.DataFrame({'Protein Name': protein_names_APR, 'Sequence_Length': sequence_length}) 
            count = 0
            while count < len(sequence_length):
                array1 = list(range(1, sequence_length[count] + 1))
                array2 = APR_residue_lists[count]
                common_elements = list(filter(lambda x: x in array2, array1))
                APR_in_total_all.append(len(common_elements)*100/sequence_length[count])
                count += 1
            protein_names_all_2params = protein_names_APR
        array1_candis = ['DFR_action', 'Sbind_action' , 'APR_action', 'MOTIF_action']
        array2_candis = ['DFR_action', 'Sbind_action' , 'APR_action', 'MOTIF_action']
        all_2params = [[], []]
        for i in array1_candis:
            if i in globals():
                i_bool = globals()[i]
            if i_bool == False:
                continue
            for j in array2_candis:
                array1_in_array2 = [[], []]
                j_flag = False
                if j in globals():
                    j_bool = globals()[j]
                if j_bool == False or i == j: continue
                count = 0
                while count < len(sequence_length):
                    if i == 'DFR_action': 
                        array1 = DFR_residue_lists[count]
                        param1 = 'DFR'
                    elif i == 'Sbind_action': 
                        array1 = Sbind_residue_lists[count]
                        param1 = 'Sbind'
                    elif i == 'APR_action': 
                        array1 = APR_residue_lists[count]
                        param1 = 'APR'
                    elif i == 'MOTIF_action': 
                        array1 = residue_lists[count]
                        param1 = "MOTIF"
                    if j == 'DFR_action': 
                        array2 = DFR_residue_lists[count]
                        param2 = "DFR"
                    elif j == 'Sbind_action': 
                        array2= Sbind_residue_lists[count]
                        param2 = "Sbind"
                    elif j == 'APR_action': 
                        array2 = APR_residue_lists[count]
                        param2 = "APR"
                    elif j == 'MOTIF_action': 
                        array2 = residue_lists[count]
                        param2 = "MOTIF"
                    common_elements = list(filter(lambda x: x in array2, array1))
                    if j_flag == False:
                        array1_in_array2[0].append(f"{param1}_in_{param2}")
                        j_flag = True
                    result = len(common_elements)*100/len(array2)
                    rounded_result = round(result, 2)
                    array1_in_array2[1].append(rounded_result)
                    count += 1
                all_2params[0].append(array1_in_array2[0])
                all_2params[1].append(array1_in_array2[1])

        data_dict = {
            'Protein Name': protein_names_all_2params,
        }
        for i in [DFR_in_total_all, Sbind_in_total_all, APR_in_total_all, MOTIF_in_total_all]:
            if len(i) > 0:
                if i == DFR_in_total_all: col_name = 'DFR_in_Sequence'
                elif i == Sbind_in_total_all: col_name = 'Sbind_in_Sequence'
                elif i == APR_in_total_all: col_name = 'APR_in_Sequence'
                elif i == MOTIF_in_total_all: col_name = 'MOTIF_in_Sequence'
                data_dict[f'{col_name}'] = i
        for i, j in zip(all_2params[0], all_2params[1]):
            col_name = ""
            for a in i:
                col_name += a
            data_dict[f'{col_name}'] = j

        all_2params_df = pd.DataFrame(data_dict)
        output_file_path = 'new_protein_data_all_2params.csv'
        all_2params_df.to_csv(output_file_path, index=False)
        seq_output_file_path = 'new_protein_data_seq.csv'
        seq_len_df.to_csv(seq_output_file_path, index=False)

        file_pattern = os.path.join(directory, 'new_protein_data*.csv')
        file_list = glob.glob(file_pattern)
        files_with_dates = [(file, os.path.getctime(file)) for file in file_list]
        files_sorted_by_date = sorted(files_with_dates, key=lambda x: x[1])
        i = 0
        switch = 0    
        if len(files_sorted_by_date) == 1:
                output_file =  directory + "\\" +  csv_plot_output + ".csv" 
                os.rename(files_sorted_by_date[0][0] , output_file)
        else:
            while i < len(files_sorted_by_date):
                if switch == 0:
                    locfile1 = files_sorted_by_date[i][0]
                    locfile2 = files_sorted_by_date[i+1][0]
                    switch = 1
                    i += 1
                elif switch == 1:
                    locfile1 = f'new_protein_data_temp_{i-1}.csv'
                    locfile2 = files_sorted_by_date[i][0]
                output_file = f'new_protein_data_temp_{i}.csv'
                if i == len(files_sorted_by_date) -1:
                    output_file =  directory + "\\" +  csv_plot_output + ".csv" 
                with open(locfile1, 'r', newline='') as file1:
                    reader1 = list(csv.reader(file1))
                with open(locfile2, 'r', newline='') as file2:
                    reader2 = list(csv.reader(file2))
                max_rows = max(len(reader1), len(reader2))
                len_row1 = len(reader1[0]) if reader1 else 0
                len_row2 = len(reader2[0]) if reader2 else 0
                while len(reader1) < max_rows:
                    reader1.append([""] * len_row1)
                while len(reader2) < max_rows:
                    reader2.append([""] * len_row2)
                with open(output_file, 'w', newline='') as out_file:
                    writer = csv.writer(out_file)
                    for row1, row2 in zip(reader1, reader2):
                        merged_row = row1 + row2
                        writer.writerow(merged_row)
                i += 1          
        files_to_delete = glob.glob(os.path.join(directory, 'new_protein_data*.csv'))
        for file_path in files_to_delete:
            try:
                os.remove(file_path)
                print(f"Deleted file: {file_path}")
            except Exception as e:
                print(f"Error deleting file {file_path}: {e}")
def plotter(residue_map, Sbind_DFR_map, APR_map, directory, csv_plot_output, DFR_action, Sbind_action, Aggregate_action, MOTIF_action, col_DFR, col_Sbind, col_MOTIF, col_APR):
        legend_handles = [Line2D([0], [0], color='black', lw=2, label='Sequence length')]
        Sbind_map = Sbind_DFR_map[0]
        DFR_map = Sbind_DFR_map[1]
        unique_weights_res = []
        unique_weights_DFR = []
        unique_weights_Sbind = []
        unique_weights_APR = []
        if residue_map != False: X_map = residue_map
        else:     X_map = {'weight_lists': [[0]]}
        if DFR_map != False: Y_map = DFR_map
        else:     Y_map = {'weight_lists': [[0]]}
        if Sbind_map != False: Z_map = Sbind_map
        else:     Z_map = {'weight_lists': [[0]]}
        if APR_map != False: W_map = APR_map
        else:     W_map = {'weight_lists': [[0]]}
        x_MAP = []
        y_MAP = []
        z_MAP = []
        w_MAP = []
        for array_of_arrays, out in zip([X_map['weight_lists'], Y_map['weight_lists'], Z_map['weight_lists'], W_map['weight_lists']], [x_MAP, y_MAP, z_MAP, w_MAP]):
            for sub_array in array_of_arrays:
                for value in sub_array:
                    out.append(value)
        for i, j, k, l in itertools.zip_longest(x_MAP, y_MAP, z_MAP, w_MAP):
            if unique_weights_res.count(i) == 0: unique_weights_res.append(i)
            if unique_weights_DFR.count(j) == 0: unique_weights_DFR.append(j)
            if unique_weights_Sbind.count(k) == 0: unique_weights_Sbind.append(k)
            if unique_weights_APR.count(l) == 0: unique_weights_APR.append(l)
        scale_determine = max(len(unique_weights_res), len(unique_weights_DFR), len(unique_weights_Sbind), len(unique_weights_APR))
        if scale_determine <= 5:
            weight_scaling_factor = 0.5
            fig, ax = plt.subplots(figsize=(12, 2))
            legend_pos = (0.5, -0.13)
        elif scale_determine <= 10:
            weight_scaling_factor = 1
            fig, ax = plt.subplots(figsize=(12, 4))
            legend_pos = (0.5, -0.1)
        elif scale_determine <=20:
            weight_scaling_factor = 2
            fig, ax = plt.subplots(figsize=(12, 6))
            legend_pos = (0.5, -0.08)
        elif scale_determine <= 30:
            weight_scaling_factor = 4
            fig, ax = plt.subplots(figsize=(12, 8))
            legend_pos = (0.5, -0.06)
        else:
            weight_scaling_factor = 7
            fig, ax = plt.subplots(figsize=(12, 10))
            legend_pos = (0.5, -0.04)
        plt.xlabel('Residue Number')
        ## PLOTS DFR and/or Sbind Residues
        pos = []
        labels = []
        if Sbind_action or DFR_action:
            if DFR_action:
                x_map = DFR_map
                legend_handles.append(Line2D([0], [0], marker='o', color='w', label='DFR Regions', markerfacecolor= col_DFR, markersize=10))
            if Sbind_action:
                x_map = Sbind_map
                legend_handles.append(Line2D([0], [0], marker='o', color='w', label='Sbind Regions', markerfacecolor= col_Sbind, markersize=10))     
            protein_names = DFR_map['protein_names']   
            Sbind_protein_name = Sbind_map['protein_names']
            length_Sbind = len(Sbind_protein_name) 
            length_DFR = len(protein_names)
            count = 0
            while count < max(length_Sbind, length_DFR) :
                end_residue = x_map['sequence_length'][count]
                if DFR_action == True and count < length_DFR:
                    DFR_y_list = []
                    for i in DFR_map['weight_lists'][count]:
                        DFR_y_list.append(weight_scaling_factor*i)    
                    DFR_x_values = DFR_map['residue_lists'][count]
                    DFR_y_values = DFR_y_list
                    plt.scatter(DFR_x_values, DFR_y_values, color= col_DFR , s = 35)
                if Sbind_action == True and count < length_Sbind:
                    Sbind_y_list = []
                    for i in Sbind_map['weight_lists'][count]:
                        Sbind_y_list.append(weight_scaling_factor*i)    
                    Sbind_x_values = Sbind_map['residue_lists'][count]
                    Sbind_y_values = Sbind_y_list
                    plt.scatter(Sbind_x_values, Sbind_y_values, color=col_Sbind, s = 20)
                temp = x_map['weight_lists'][count][0]
                pos.append(weight_scaling_factor*temp)
                protein_name = x_map['protein_names'][count]
                labels.append(protein_name)
                start_residue = 1
                #plt.hlines(y=weight_scaling_factor*temp, xmin=start_residue, xmax=end_residue, color='grey', linestyle='-', linewidth=0.3)
                count += 1
        ## PLOTS MAPPED Residues
        if MOTIF_action == True:
            len_residue_map = len(residue_map['sequence_length'])
            res_count = 0
            while res_count < len_residue_map:
                end_residue = residue_map['sequence_length'][res_count]
                y_list = []
                for i in residue_map['weight_lists'][res_count]:
                    y_list.append(weight_scaling_factor*i)
                y_values = y_list
                x_values = residue_map['residue_lists'][res_count]
                plt.scatter(x_values, y_values, color= col_MOTIF , s = 9)
                temp = residue_map['weight_lists'][res_count][0]
                pos.append(weight_scaling_factor*temp)
                protein_name = residue_map['protein_names'][res_count]
                labels.append(protein_name)
                start_residue = 1
                #plt.hlines(y=weight_scaling_factor*temp, xmin=start_residue, xmax=end_residue, color='grey', linestyle='-', linewidth=0.3)
                res_count += 1
            MOTIF_name = residue_map['MOTIF type']
            temp_name = ''
            for i in MOTIF_name[0]:
                temp_name += i
            if MOTIF_headings.count(temp_name) == 0:
                temp_name = temp_name + " Repeats"
            legend_handles.append(Line2D([0], [0], marker='o', color='w', label= temp_name , markerfacecolor= col_MOTIF, markersize=10))
        if Aggregate_action == True:
            len_APR_map = len(APR_map['protein_names'])
            APR_count = 0
            while APR_count < len_APR_map:
                end_residue = APR_map['sequence_length'][APR_count]
                y_list = []
                for i in APR_map['weight_lists'][APR_count]:
                    y_list.append(weight_scaling_factor*i)
                y_values = y_list
                x_values = APR_map['residue_lists'][APR_count]
                plt.scatter(x_values, y_values, color=col_APR, s = 4)
                temp = APR_map['weight_lists'][APR_count][0]
                pos.append(weight_scaling_factor*temp)
                protein_name = APR_map['protein_names'][APR_count]
                labels.append(protein_name)
                start_residue = 1
                plt.hlines(y=weight_scaling_factor*temp, xmin=start_residue, xmax=end_residue, color='grey', linestyle='-', linewidth=0.3)
                APR_count += 1
            legend_handles.append(Line2D([0], [0], marker='o', color='w', label= 'APR' , markerfacecolor=col_APR, markersize=10))
        
        plt.hlines(y=weight_scaling_factor*temp, xmin=start_residue, xmax=end_residue, color='grey', linestyle='-', linewidth=0.3)
        plt.hlines(y=weight_scaling_factor*temp, xmin=start_residue, xmax=end_residue, color='grey', linestyle='-', linewidth=0.3)
        plt.yticks(pos, labels, fontsize = 10, color = 'grey')
        plt.legend(handles=legend_handles, loc='upper center', bbox_to_anchor= legend_pos , fancybox=True, shadow=True, ncol=3)
        plt.tight_layout()
        output_file = directory + "\\" +  csv_plot_output + '.png'
        plt.savefig(output_file)
        plt.show()
'''
if len(sys.argv) > 1: 
    directory = sys.argv[1]
else:
    directory = input("Enter working directory: ")
if len(sys.argv) > 2: 
    protein_sequence = sys.argv[2]
else:
    protein_sequence = input("Enter file containing sequence data: ")
if len(sys.argv) > 3: 
    Sbind_DFR_sequence = sys.argv[3]
else:
    DFR_sequence = input("Enter file containing DFR data: ") 
    '''
cyan = '#00FFFF'
good_green = '#2ADE2A'
blue = '#2020d4'
red = '#de2a2a'
col_DFR = cyan
col_Sbind = red
col_MOTIF = 'blue'
col_APR = 'black'
value_cuttoff = 21.4
value_repeat_num = 5
sbind_value_cuttoff = 2.2
sbind_value_repeat_num = 10
pDP_value_cuttoff = 0.6
flag = True
while flag == True:         
    if len(sys.argv) > 4: 
        csv_plot_output = sys.argv[4]
    else:
        csv_plot_output = input("Enter output file name (9999 to terminate): ") 
    if csv_plot_output == "9999":
        break
    if len(sys.argv) > 5: 
        actions = sys.argv[5]
    else:
        actions = input("Possible actions(enter all applicable continuosly):\n1.DFR_map\n2.Sbind_map\n3.Aggregate_map\n4.Motif_map\n") 
    Sbind_action = DFR_action = MOTIF_action = APR_action = False
    if actions.count("1"):
        DFR_action = True
    if actions.count("2"):
        Sbind_action = True
    if actions.count("3"):
        APR_action = True
    if actions.count("4"):
        MOTIF_action = True
    
    aggreagate_data = r"HOX_APR"
    seq_filename = r"C:\Users\KANT4\Desktop\New folder\React\Ongoing\scatter_plot\src\xls datasets\protein_sequence.txt"
    Sbind_DFR_sequence = r"DFR_full"
    directory = r"C:\Users\KANT4\Desktop\New folder\React\Ongoing\scatter_plot\src\xls datasets\csv generated"
        
    # seq_filename = directory + "\\" +  protein_sequence + ".txt"
    Sbind_DFR_filename = directory + "\\" +  Sbind_DFR_sequence + ".txt" 
    APR_filename =  directory + "\\" +  aggreagate_data + ".txt" 
    if Sbind_action == True or DFR_action == True:
        Sbind_DFR_map = Sbind_DFR_mapper(seq_filename, Sbind_DFR_filename, sbind_value_cuttoff, sbind_value_repeat_num, pDP_value_cuttoff)
    else:
        Sbind_DFR_map = [False, False]
    if MOTIF_action == True:
        residue_map = motif_mapper(seq_filename)
    else:
        residue_map = False
    if APR_action == True:
        APR_map = APR_mapper(APR_filename, value_cuttoff, value_repeat_num)
    else:
        APR_map = False
    make_CSV(residue_map, Sbind_DFR_map, APR_map, directory, csv_plot_output, DFR_action, Sbind_action, APR_action, MOTIF_action)         
    plotter(residue_map, Sbind_DFR_map, APR_map, directory, csv_plot_output, DFR_action, Sbind_action, APR_action, MOTIF_action, col_DFR, col_Sbind, col_MOTIF, col_APR)