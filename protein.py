import pandas as pd
import re
from sqlalchemy import create_engine, Column, Integer, String, Float, MetaData, Table
from app import engine

# Read data from the file into a Pandas DataFrame
file_path = 'INDEX_core_name.2013'
data_list = []
pattern = re.compile(r'\((.*?)\)')

with open(file_path, 'r') as file:
    lines = [line.strip() for line in file if not line.startswith('#')]
    for line in lines:
        fields = line.split('$')
        pdb_code, release_year, uniprot_id, protein_name = fields
        
        data_list.append({
            'pdb_code': pdb_code,
            'release_year': release_year,
            'uniprot_id': uniprot_id,
            'protein_name': protein_name,
        })


df = pd.DataFrame(data_list)

df.to_sql('proteins', con=engine, if_exists='replace', index=False)