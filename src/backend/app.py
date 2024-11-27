from flask import Flask, request, jsonify, make_response
from flask_cors import CORS
import csv
import re
import io
import logging

app = Flask(__name__)
CORS(app)

# Setup logging
logging.basicConfig(level=logging.DEBUG)

# In-memory data store
data_store = {}

# File paths (update these as needed)
seq_filename = r"..\xls datasets\protein_sequence.txt"
ranks_scores_file_path = r"..\xls datasets\rank_score.csv"


rank_score_int_flag = False  # Update flag as needed

def read_columns_to_arrays(file_path, int_flag):
    first_column_values = []
    second_column_values = []
    first_column_values_as_string = []
    with open(file_path, mode='r', newline='', encoding='utf-8-sig') as file:
        reader = csv.reader(file)
        flag = True
        for row in reader:
            if flag:
                flag = False
                continue
            if not int_flag:
                first_column_values.append(row[0])
            else:
                first_column_values.append(int(row[0]))
            first_column_values_as_string.append(f"{row[0]}")
            second_column_values.append(rf'{row[1]}')
    return first_column_values, second_column_values, first_column_values_as_string

def motif_mapper(seq_filename, Ranked_Weights, Scores, input_value):
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
        data['ranked_weight_lists'].append(ranked_weights)
        data['score_lists'].append(scores)
        data['MOTIF type'].append(MOTIF_type)
        data['sequence_length'].append(seq_len)

    motif_map = create_residue_map()
    MOTIF_type = input_value  # Now using the input_value passed as a parameter
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
                    residue_list.append(index)
                    weight_list.append(temp)
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

    def save_motif_map_to_csv(motif_map):
        headers = ['Protein Name', 'Residues', 'Amino Acids', 'Weights', 'Ranked Weights', 'Scores', 'Sequence Length']
        rows = zip(
            motif_map['protein_names'],
            motif_map['residue_lists'],
            motif_map['amino_acid_lists'],
            motif_map['weight_lists'],
            motif_map['ranked_weight_lists'],
            motif_map['score_lists'],
            motif_map['sequence_length']
        )
        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerow(headers)
        for row in rows:
            writer.writerow(row)
        return output.getvalue()

    csv_content = save_motif_map_to_csv(motif_map)
    return csv_content

@app.route('/generate-data', methods=['POST'])
def generate_data():
    input_value = request.json.get('input')
    if not input_value:
        return jsonify({'error': 'No input provided'}), 400

    label = f"{input_value}"
    logging.debug(f"Generated label: {label}")
    Ranked_Weights, Scores, _ = read_columns_to_arrays(ranks_scores_file_path, rank_score_int_flag)
    csv_content = motif_mapper(seq_filename, Ranked_Weights, Scores, input_value)
    logging.debug(f"Generated CSV content for label: {label}")

    data_store[label] = csv_content
    return jsonify({'label': label})

@app.route('/labels', methods=['GET'])
def get_labels():
    try:
        labels = list(data_store.keys())
        return jsonify({'labels': labels})
    except Exception as e:
        logging.error(f"Error fetching labels: {str(e)}")
        return jsonify({'error': str(e)}), 500

@app.route('/data/<label>', methods=['GET'])
def get_data(label):
    label = label.replace('.csv', '')  # Ensure consistency
    if label in data_store:
        csv_content = data_store[label]
        response = make_response(csv_content)
        response.headers['Content-Disposition'] = f'attachment; filename={label}.csv'
        response.headers['Content-Type'] = 'text/csv'
        return response
    else:
        logging.error(f"Label not found: {label}")
        return jsonify({'error': 'Label not found'}), 404

@app.route('/download-csv/<label>')
def download_csv(label):
    label = label.replace('.csv', '')  # Ensure consistency
    if label in data_store:
        csv_content = data_store[label]
        response = make_response(csv_content)
        response.headers['Content-Disposition'] = f'attachment; filename={label}.csv'
        response.headers['Content-Type'] = 'text/csv'
        return response
    else:
        logging.error(f"Label not found: {label}")
        return jsonify({'error': 'Label not found'}), 404

if __name__ == '__main__':
    app.run(debug=True)
