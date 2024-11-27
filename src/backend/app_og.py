from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
import csv
import os
import uuid

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

DATA_DIR = 'data'

if not os.path.exists(DATA_DIR):
    os.makedirs(DATA_DIR)

@app.route('/generate-data', methods=['POST'])
def generate_data():
    input_value = request.json.get('input')
    if not input_value:
        return jsonify({'error': 'No input provided'}), 400

    label = f"{input_value}_{uuid.uuid4()}"
    data = generate_random_data(label)

    file_path = os.path.join(DATA_DIR, f"{label}.csv")
    save_data_to_csv(file_path, data)

    return jsonify({'label': label})

@app.route('/data/<filename>')
def get_data(filename):
    file_path = os.path.join(DATA_DIR, filename)
    if os.path.exists(file_path):
        return send_file(file_path, as_attachment=True)
    else:
        return jsonify({'error': 'File not found'}), 404

def generate_random_data(label):
    import random
    data = [{'x': random.uniform(-100, 100), 'y': random.uniform(-100, 100), 'id': i} for i in range(100)]
    return data

def save_data_to_csv(file_path, data):
    with open(file_path, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=['x', 'y', 'id'])
        writer.writeheader()
        for row in data:
            writer.writerow(row)

if __name__ == '__main__':
    app.run(debug=True) 