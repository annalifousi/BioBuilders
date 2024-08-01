from flask import Flask, request, render_template, jsonify
import sqlite3
import pandas as pd
from sqlalchemy import create_engine

app = Flask(__name__)

# Define database path
db_path = 'ner_results.db'

# Initialize the database
def initialize_db():
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS ner_entities (
        Entity_EDC TEXT,
        Entity_TARGET TEXT,
        Count INTEGER,
        articles INTEGER
    )
    ''')
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

        # Connect to the database and filter based on uploaded file
        engine = create_engine(f'sqlite:///{db_path}')
        query = f"""
        SELECT * FROM sentences 
        WHERE endocrine_disrupting_chemical IN ({','.join(['?']*len(edcs_or_targets))}) 
        OR target IN ({','.join(['?']*len(edcs_or_targets))})
        """
        conn = sqlite3.connect(db_path)
        df = pd.read_sql_query(query, conn, params=edcs_or_targets + edcs_or_targets)
        conn.close()

        # Convert the result to a dictionary for JSON response
        result = df.to_dict(orient='records')

        return jsonify(result)
    
    return jsonify({'error': 'File upload failed'}), 400

if __name__ == '__main__':
    app.run(debug=True)
