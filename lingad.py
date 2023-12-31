import pandas as pd
import re
from sqlalchemy import create_engine, Column, Integer, String, Float, MetaData, Table
from app import engine

# Read data from the file into a Pandas DataFrame
file_path = 'INDEX_core_data.2013'
data_list = []
pattern = re.compile(r'\((.*?)\)')

with open(file_path, 'r') as file:
    lines = [line.strip() for line in file if not line.startswith('#')]
    for line in lines:
        fields = line.split('$')
        print(fields)
        pdb_code, resolution, release_year,log, binding_data, reference, ligand_name = fields
        math = pattern.search(ligand_name)
        ligand_name = math.group(1) if math else None
        data_list.append({
            'pdb_code': pdb_code,
            'resolution': resolution,
            'release_year': release_year,
            'logKd/Ki' : log,
            'binding_data': binding_data,
            'reference': reference,
            'ligand_name': ligand_name
        })


df = pd.DataFrame(data_list)

df.to_sql('lingads', con=engine, if_exists='replace', index=False)