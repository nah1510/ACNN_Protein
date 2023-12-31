import os
import numpy as np
import pandas as pd

import tempfile

from rdkit import Chem
from rdkit.Chem import AllChem
import deepchem as dc

from deepchem.utils import download_url, load_from_disk 

data_dir = dc.utils.get_data_dir()
dataset_file = os.path.join(data_dir, "pdbbind_core_df.csv.gz")

if not os.path.exists(dataset_file):
    print('File does not exist. Downloading file...')
    download_url("https://s3-us-west-1.amazonaws.com/deepchem.io/datasets/pdbbind_core_df.csv.gz")
    print('File downloaded...')
    
    
raw_dataset = load_from_disk(dataset_file)
raw_dataset = raw_dataset[['pdb_id', 'smiles', 'label']]

from openmm.app import PDBFile
from pdbfixer import PDBFixer
import mdtraj as md
import nglview
from deepchem.utils.vina_utils import prepare_inputs
from IPython.display import display, Image

for i in range(171, len(raw_dataset)): 
    pdbid = raw_dataset['pdb_id'].iloc[i]
    ligand = raw_dataset['smiles'].iloc[i]
    
    print(i)
    print(pdbid)
    
    # some pdb code error
    if pdbid in ["2zjw","1e66","1r5y","3su2","1hfs"] :
        continue
    fixer = PDBFixer(pdbid=pdbid)
    PDBFile.writeFile(fixer.topology, fixer.positions, open('templates/display/%s.pdb' % (pdbid), 'w'))
    p, m = None, None
    # fix protein, optimize ligand geometry, and sanitize molecules
    try:
        print(pdbid, ligand)
        
        p, m = prepare_inputs('templates/display/%s.pdb' % (pdbid), ligand)
    except:
        print('%s failed PDB fixing' % (pdbid)) 
    if p and m:  # protein and molecule are readable by RDKit
        Chem.rdmolfiles.MolToPDBFile(p, 'templates/display/protein_%s.pdb' % (pdbid))
        Chem.rdmolfiles.MolToPDBFile(m, 'templates/display/ligand_%s.pdb' % (pdbid))
    
    
        
    protein_mdtraj = md.load_pdb( 'templates/display/protein_%s.pdb' % (pdbid))
    ligand_mdtraj = md.load_pdb('templates/display/ligand_%s.pdb' % (pdbid))
    p = nglview.show_mdtraj(protein_mdtraj)
    l = nglview.show_mdtraj(ligand_mdtraj)
    nglview.write_html("templates/display/protein_%s.html" % (pdbid),[p])
    nglview.write_html("templates/display/ligand_%s.html" % (pdbid),[l])
    