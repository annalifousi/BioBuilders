from flask import Flask, render_template
from sqlalchemy import create_engine, Column, Integer, String
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.exc import OperationalError

app = Flask(__name__)

# Database setup
DATABASE_URI = 'sqlite:///ner_results.db'
engine = create_engine(DATABASE_URI, echo=True)  # Enable SQL echo for debugging
Session = sessionmaker(bind=engine)
session = Session()
Base = declarative_base()

# Define the table structure using SQLAlchemy ORM
class Sentences(Base):
    __tablename__ = 'sentences'
    activity = Column('activity', String, primary_key=True)
    endocrine_disrupting_chemical = Column('endocrine_disrupting_chemical', String, primary_key=True)
    target = Column('target', String, primary_key=True)
    counts = Column('counts', Integer)
    articles = Column('articles', Integer)

@app.route('/')
def index():
    try:
        # Fetch data from the database
        results = session.query(Sentences).all()
        return render_template('index.html', results=results)
    except OperationalError as e:
        return f"An error occurred: {e}"

if __name__ == '__main__':
    app.run(debug=True)
