from flask import Flask, request, render_template, jsonify
from flask_cors import CORS
import sqlite3
import pandas as pd
from sqlalchemy import create_engine

app = Flask(__name__)
CORS(app)  # Enable CORS

# Define database path
db_path = 'ner_results.db'

# Initialize the database
def initialize_db():
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Create or update the ner_entities table with the correct schema
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS ner_entities (
        endocrine_disrupting_chemical TEXT,
        activity TEXT,
        target TEXT,
        counts INTEGER,
        articles INTEGER
    )
    ''')
    
    # Create or update the sentences table with the correct schema
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS sentences (
        activity TEXT,
        endocrine_disrupting_chemical TEXT,
        target TEXT,
        counts INTEGER,
        articles INTEGER,
        PRIMARY KEY (activity, endocrine_disrupting_chemical, target)
    )
    ''')
    
    conn.commit()
    conn.close()

# Initialize the database
initialize_db()

# Route for the home page with upload form
@app.route('/')
def index():
    return render_template('index.html')

# Route to handle search queries
@app.route('/search', methods=['GET'])
def search():
    query = request.args.get('query', '')

    conn = sqlite3.connect(db_path)
    query_str = """
    SELECT * FROM sentences 
    WHERE endocrine_disrupting_chemical LIKE ? 
    OR target LIKE ?
    ORDER BY endocrine_disrupting_chemical, activity, target
    """
    params = [f"%{query}%", f"%{query}%"]
    df = pd.read_sql_query(query_str, conn, params=params)
    conn.close()

    return jsonify(df.to_dict(orient='records'))

# Route to handle file upload and filtering
@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'}), 400
    
    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'No selected file'}), 400
    
    if file:
        # Read the uploaded file
        edcs_or_targets = file.read().decode('utf-8').splitlines()
        
        if not edcs_or_targets:
            return jsonify({'error': 'Empty file content'}), 400
        
        # Connect to the database and filter based on uploaded file
        query_placeholders = ','.join(['?'] * len(edcs_or_targets))
        query = f"""
        SELECT * FROM sentences 
        WHERE endocrine_disrupting_chemical IN ({query_placeholders}) 
        OR target IN ({query_placeholders})
        ORDER BY endocrine_disrupting_chemical, activity, target
        """
        conn = sqlite3.connect(db_path)
        params = edcs_or_targets * 2  # List should be twice the length for both placeholders
        df = pd.read_sql_query(query, conn, params=params)
        conn.close()

        # Convert the result to a dictionary for JSON response
        result = df.to_dict(orient='records')

        return jsonify(result)
    
    return jsonify({'error': 'File upload failed'}), 400

# Route to provide data for network visualization
@app.route('/network', methods=['GET'])
def network():
    target = request.args.get('target', '')

    conn = sqlite3.connect(db_path)
    query_str = """
    SELECT * FROM sentences 
    WHERE target LIKE ?
    ORDER BY endocrine_disrupting_chemical, activity, target
    """
    params = [f"%{target}%"]
    df = pd.read_sql_query(query_str, conn, params=params)
    conn.close()

    nodes = [{'data': {'id': target, 'label': target, 'size': 50}}]  # Fixed size for the central node
    edges = []

    for _, row in df.iterrows():
        edc = row['endocrine_disrupting_chemical']
        counts = row['counts']
        nodes.append({'data': {'id': edc, 'label': edc, 'size': counts}})
        edges.append({'data': {'source': target, 'target': edc}})

    return jsonify({'nodes': nodes, 'edges': edges})

if __name__ == '__main__':
    app.run(debug=True)
