import deepchem as dc
import os

import numpy as np
# import tensorflow as tf

import matplotlib.pyplot as plt
import joblib

from rdkit import Chem

from load_input import load_input
from deepchem.models import AtomicConvModel
from deepchem.feat import AtomicConvFeaturizer

f1_num_atoms = 100  # maximum number of atoms to consider in the ligand
f2_num_atoms = 1000  # maximum number of atoms to consider in the protein
max_num_neighbors = 12  # maximum number of spatial neighbors for an atom

acf = AtomicConvFeaturizer(frag1_num_atoms=f1_num_atoms,
                      frag2_num_atoms=f2_num_atoms,
                      complex_num_atoms=f1_num_atoms+f2_num_atoms,
                      max_num_neighbors=max_num_neighbors,
                      neighbor_cutoff=4)

tasks, [train, val, input], transformers = load_input(featurizer=acf,
                                             save_dir='.',
                                             data_dir='.',
                                             pocket=True,
                                             reload=False,
                                             protein_pdbcode="4des",
                                             ligand_pdbcode="4djv",
                                             set_name='core')

print(input)
