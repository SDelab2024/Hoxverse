import pandas as pd
import csv
import re
import os

def read_columns_to_arrays(file_path, int_flag):
    first_column_values = []
    second_column_values = []
    first_column_values_as_string = []
    with open(file_path, mode='r', newline='', encoding='utf-8-sig') as file:
        reader = csv.reader(file)
        flag = True
        for row in reader:
            if flag:  # skips line with heading
                flag = False
                continue
            if int_flag:
                first_column_values.append(int(row[0]))
            else:
                first_column_values.append(row[0])
            first_column_values_as_string.append(f"{row[0]}")
            second_column_values.append(rf'{row[1]}')

    return first_column_values, second_column_values, first_column_values_as_string

def create_residue_map():
    return {
        'protein_names': [],
        'residue_lists': [],
        'amino_acid_lists': [],
        'weight_lists': [],
        'ranked_weight_lists': [],
        'score_lists': [],
        'MOTIF type': [],
        'sequence_length': []
    }

def add_residue_map(data, protein_name, residues, amino_acids, weights, ranked_weights, scores, MOTIF_type, seq_len):
    data['protein_names'].append(protein_name)
    data['residue_lists'].append(residues)
    data['amino_acid_lists'].append(amino_acids)
    data['weight_lists'].append(weights)
    data['ranked_weight_lists'].append(ranked_weights)
    data['score_lists'].append(scores)
    data['MOTIF type'].append(MOTIF_type)
    data['sequence_length'].append(seq_len)

def motif_mapper(seq_filename, Ranked_Weights, Scores):
    MOTIF_ranked_weights = Ranked_Weights
    MOTIF_Scores = Scores
    
    motif_map = create_residue_map()
    motif_user_input = input("Enter custom motif: ")
    MOTIF_type = motif_user_input
    motif_map['MOTIF type'].append(MOTIF_type)
    with open(seq_filename, 'r') as file:
        lines = file.readlines()

    num_seqs = sum(1 for line in lines if line.startswith('>'))

    temp = num_seqs
    rank_count = 0
    for line in lines:
        if line.startswith('>'):
            no_matches_flag = True
            index = 0
            protein_name = line.split('|')[2].split('_')[0]
            residue_list = []
            weight_list = []
            ranked_weight_list = []
            score_list = []
            amino_acid_list = []
            for next_line in lines[lines.index(line) + 1:]:
                if next_line.startswith('>'):
                    break
                for char in next_line.strip():
                    index += 1
                    amino_acid_list.append(char)
                    residue_list.append(index)  # Assuming index represents residue number
                    weight_list.append(temp)   # Assuming temp represents some weight
                    ranked_weight_list.append(MOTIF_ranked_weights[rank_count])
                    score_list.append(MOTIF_Scores[rank_count])
            sequence = ''.join(amino_acid_list)
            seq_len = len(sequence)
            matches = re.finditer(MOTIF_type, sequence)
            for match in matches:
                no_matches_flag = False
                start_index = match.start()
                end_index = match.end()
                motif_map['protein_names'].append(protein_name)
                motif_map['residue_lists'].append(residue_list[start_index:end_index])
                motif_map['amino_acid_lists'].append(amino_acid_list[start_index:end_index])
                motif_map['weight_lists'].append(weight_list[start_index:end_index])
                motif_map['sequence_length'].append(seq_len)
                motif_map['ranked_weight_lists'].append(ranked_weight_list[start_index:end_index])
                motif_map['score_lists'].append(score_list[start_index:end_index])
            if no_matches_flag:
                motif_map['protein_names'].append(protein_name)
                motif_map['residue_lists'].append([0])
                motif_map['amino_acid_lists'].append([" "])
                motif_map['weight_lists'].append([temp])
                motif_map['sequence_length'].append(seq_len)
                motif_map['ranked_weight_lists'].append([MOTIF_ranked_weights[rank_count]])
                motif_map['score_lists'].append([MOTIF_Scores[rank_count]])
            temp -= 1
            rank_count += 1

    # Removing duplicates and merging data
    names = motif_map['protein_names']
    all_res_lists = motif_map['residue_lists']
    all_aa_lists = motif_map['amino_acid_lists']
    all_weight_lists = motif_map['weight_lists']
    all_ranked_weight_lists = motif_map['ranked_weight_lists']
    all_score_lists = motif_map['score_lists']
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
        ranked_weight_list = []
        score_list = []
        seq = []
        while name == prev_name:
            for a, b, c, d, e in zip(all_res_lists[index], all_aa_lists[index], all_weight_lists[index], all_ranked_weight_lists[index], all_score_lists[index]):
                res_list.append(a)
                aa_list.append(b)
                weight_list.append(c)
                ranked_weight_list.append(d)
                score_list.append(e)
            if index < len(names) - 1:
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
            del ranked_weight_list[pindex]
            del score_list[pindex]
        motif_map['protein_names'].append(curr_name)
        motif_map['residue_lists'].append(res_list)
        motif_map['amino_acid_lists'].append(aa_list)
        motif_map['weight_lists'].append(weight_list)
        motif_map['sequence_length'].append(seq)
        motif_map['ranked_weight_lists'].append(ranked_weight_list)
        motif_map['score_lists'].append(score_list)
        prev_name = name

    return motif_map

def make_CSV(residue_map, directory, csv_plot_output, MOTIF_action):
    MOTIF_in_total_all = []
    if MOTIF_action:
        protein_names = residue_map['protein_names']
        residue_lists = residue_map['residue_lists']
        amino_acid_lists = residue_map['amino_acid_lists']
        weight_lists = residue_map['weight_lists']
        ranked_weight_lists = residue_map['ranked_weight_lists']
        score_lists = residue_map['score_lists']
        data = {'Protein Name': [], 'Residues': [], 'Amino_acids': [], 'Weights': [], 'Ranked_Weights': [], 'Scores': []}
        for protein_name, residues, amino_acids, weights, ranked_weights, scores in zip(protein_names, residue_lists, amino_acid_lists, weight_lists, ranked_weight_lists, score_lists):
            data['Protein Name'].extend([protein_name] * len(residues))
            data['Residues'].extend(residues)
            data['Amino_acids'].extend(amino_acids)
            data['Weights'].extend(weights)
            data['Ranked_Weights'].extend(ranked_weights)
            data['Scores'].extend(scores)
        df = pd.DataFrame(data)
        output_file_path = os.path.join(directory, 'new_protein_data.csv')
        df.to_csv(output_file_path, index=False)
        
        MOTIF_name = residue_map['MOTIF type']
        temp_name = ''.join(MOTIF_name)
        if temp_name in MOTIF_headings:
            j = MOTIF_headings.index(temp_name)
            MOTIF_definition = common_MOTIFs[j]
        else:
            MOTIF_definition = MOTIF_name
            MOTIF_name = temp_name + " Repeats"
        MOTIF_df = pd.DataFrame.from_dict(residue_map, orient='index').transpose()
        seq_len_df = MOTIF_df.pop('sequence_length')

        for count in range(len(protein_names)):
            sequence = ''.join(amino_acid_lists[count])
            motif = MOTIF_definition
            matches = re.finditer(motif, sequence)
            array1 = [match.start() for match in matches]
            array2 = residue_lists[count]
            common_elements = list(set(array1).intersection(array2))
            temp = len(common_elements) / seq_len_df[count] * 100
            MOTIF_in_total_all.append(temp)

        new_df = pd.concat([df, seq_len_df, pd.DataFrame({'MOTIF Count / Seq_Length (%)': MOTIF_in_total_all})], axis=1)
        output_file_path = os.path.join(directory, csv_plot_output)
        new_df.to_csv(output_file_path, index=False)
    else:
        output_file_path = os.path.join(directory, 'new_protein_data.csv')
        df = pd.DataFrame.from_dict(residue_map, orient='index').transpose()
        df.to_csv(output_file_path, index=False)

if __name__ == "__main__":
    directory = r"C:\Users\KANT4\Desktop\New folder\React\Ongoing\scatter_plot\src\xls datasets"
    protein_sequence = "protein_sequence.txt"
    ranks_scores_file_path = os.path.join(directory, "rank_score.csv")
    seq_filename = os.path.join(directory, protein_sequence)

    # seq_filename = input("Enter sequence filename: ")
    # Ranked_Weights = [float(x) for x in input("Enter ranked weights separated by space: ").split()]
    # Scores = [float(x) for x in input("Enter scores separated by space: ").split()]
    rank_score_int_flag = True

    directory = input("Enter directory to save files: ")
    csv_plot_output = input("Enter CSV output filename: ")
    MOTIF_action = input("Include MOTIF calculation? (yes/no): ").strip().lower() == 'yes'
    Ranked_Weights, Scores, tempsrsg = read_columns_to_arrays(ranks_scores_file_path, rank_score_int_flag)

    residue_map = motif_mapper(seq_filename, Ranked_Weights, Scores)
    make_CSV(residue_map, directory, csv_plot_output, MOTIF_action)
