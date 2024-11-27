import pandas as pd
import glob
import os
import csv
import re
def read_columns_to_arrays(file_path, int_flag):
    first_column_values = []
    second_column_values = []
    first_column_values_as_string = []
    with open(file_path, mode='r', newline='', encoding='utf-8-sig') as file:
        reader = csv.reader(file)
        # Iterate over each row in the CSV 
        flag = True
        for row in reader:
            if flag == True: ##skips line with heading
                flag = False
                continue
            if int_flag == False:  first_column_values.append(row[0])
            else: first_column_values.append(int(row[0]))
            first_column_values_as_string.append(f"{row[0]}")
            second_column_values.append(rf'{row[1]}')
            
    return first_column_values, second_column_values, first_column_values_as_string
common_MOTIFs = []
MOTIF_headings = []
def motif_mapper(seq_filename, Ranked_Weights, Scores):
    MOTIF_ranked_weights = Ranked_Weights
    MOTIF_Scores = Scores
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
        data['weight_lists'].append(weights)
        data['ranked_weight_lists'].append(ranked_weights)
        data['score_lists'].append(scores)
        data['MOTIF type'].append(MOTIF_type)
        data['sequence_length'].append(seq_len)

    motif_map = create_residue_map()  # Assuming you have a function create_residue_map() defined
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
            if no_matches_flag == True:
                motif_map['protein_names'].append(protein_name)
                motif_map['residue_lists'].append([0])
                motif_map['amino_acid_lists'].append([" "])
                motif_map['weight_lists'].append([temp])
                motif_map['sequence_length'].append(seq_len)
                motif_map['ranked_weight_lists'].append([MOTIF_ranked_weights[rank_count]])
                motif_map['score_lists'].append([MOTIF_Scores[rank_count]])
            temp -= 1
            rank_count += 1
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
    print(motif_map)
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
        data = {
            'Protein Name': [], 
            'Residues': [], 
            'Amino_acids': [], 
            'Weights': [], 
            'Ranked_Weights': [], 
            'Scores': []
        }
        for protein_name, residues, amino_acids, weights, ranked_weights, scores in zip(
            protein_names, residue_lists, amino_acid_lists, weight_lists, ranked_weight_lists, score_lists
        ):
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
        temp_name = ''.join(MOTIF_name[0])
        if temp_name in MOTIF_headings:
            j = MOTIF_headings.index(temp_name)
            MOTIF_definition = common_MOTIFs[j]
        else:
            MOTIF_definition = MOTIF_name
            temp_name += " Repeats"
        MOTIF_df = pd.DataFrame(MOTIF_definition, columns=[temp_name])
        output_file_path = os.path.join(directory, 'new_protein_data_MOTIF.csv')
        MOTIF_df.to_csv(output_file_path, index=False)
        
        sequence_length = residue_map['sequence_length']
        temp_seq_len = [b for a in sequence_length for b in a]
        sequence_length = temp_seq_len
        seq_len_df = pd.DataFrame({
            'Protein Name': protein_names, 
            'Sequence_Length': sequence_length
        })
        seq_output_file_path = os.path.join(directory, 'new_protein_data_seq.csv')
        seq_len_df.to_csv(seq_output_file_path, index=False)

        count = 0
        while count < len(sequence_length):
            array1 = list(range(1, sequence_length[count] + 1))
            array2 = residue_lists[count]
            common_elements = list(filter(lambda x: x in array2, array1))
            MOTIF_in_total_all.append(len(common_elements) * 100 / sequence_length[count])
            count += 1
        
    file_pattern = os.path.join(directory, 'new_protein_data*.csv')
    file_list = glob.glob(file_pattern)
    files_with_dates = [(file, os.path.getctime(file)) for file in file_list]
    files_sorted_by_date = sorted(files_with_dates, key=lambda x: x[1])
    
    i = 0
    switch = 0    
    if len(files_sorted_by_date) == 1:
        output_file = os.path.join(directory, csv_plot_output + ".csv") 
        os.rename(files_sorted_by_date[0][0], output_file)
    else:
        while i < len(files_sorted_by_date):
            if switch == 0:
                locfile1 = files_sorted_by_date[i][0]
                locfile2 = files_sorted_by_date[i + 1][0]
                switch = 1
                i += 1
            elif switch == 1:
                locfile1 = f'new_protein_data_temp_{i-1}.csv'
                locfile2 = files_sorted_by_date[i][0]
            
            output_file = f'new_protein_data_temp_{i}.csv'
            if i == len(files_sorted_by_date) - 1:
                output_file = os.path.join(directory, csv_plot_output + ".csv")
                
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

value_cuttoff = 21.4
value_repeat_num = 5
sbind_value_cuttoff = 2.2
sbind_value_repeat_num = 10
pDP_value_cuttoff = 0.6
       
directory = r"C:\Users\KANT4\Desktop\New folder\React\Ongoing\scatter_plot\src\xls datasets"
protein_sequence = "protein_sequence.txt"
ranks_scores_file_path = os.path.join(directory, "rank_score.csv")

rank_score_int_flag = True

# Assuming you have the necessary function definitions before this code
Ranked_Weights, Scores, tempsrsg = read_columns_to_arrays(ranks_scores_file_path, rank_score_int_flag)

seq_filename = os.path.join(directory, protein_sequence)
residue_map = motif_mapper(seq_filename, Ranked_Weights, Scores)

csv_plot_output = "Protein"
MOTIF_action = True
make_CSV(residue_map, directory, csv_plot_output, MOTIF_action)       
