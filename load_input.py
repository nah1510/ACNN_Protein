"""
PDBBind dataset loader.
"""
import os
import numpy as np
import pandas as pd

import deepchem as dc
from loader import TransformerGenerator, _MolnetLoader
from deepchem.data import Dataset
from typing import List, Optional, Tuple, Union

DATASETS_URL = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/"
PDBBIND_URL = DATASETS_URL + "pdbbindv2019/"
PDBBIND_TASKS = ['-logKd/Ki']


class _PDBBindInputLoader(_MolnetLoader):

  def __init__(self,
               *args,
               pocket: bool = True,
               set_name: str = 'core',
               **kwargs):
    super(_PDBBindInputLoader, self).__init__(*args, **kwargs)
    self.pocket = pocket
    self.set_name = set_name
    if set_name == 'general':
      self.name = 'pdbbind_v2019_other_PL'  # 'general' set folder name
    elif set_name == 'refined':
      self.name = 'pdbbind_v2019_refined'
    elif set_name == 'core':
      self.name = 'pdbbind_v2013_core_set'

  def create_dataset(self,protein_pdbcode: str,ligand_pdbcode: str) -> Dataset:
    # get pdb and sdf filenames, labels and pdbids

    labels =[ 0,0,0,0,0,0]
    pdbs = ['x','x','x','x','x','x']
    protein_files = ['./v2013-core/%s/%s_pocket.pdb' % (protein_pdbcode,protein_pdbcode),
                     './v2013-core/%s/%s_pocket.pdb' % (protein_pdbcode,protein_pdbcode),
                     './v2013-core/%s/%s_pocket.pdb' % (protein_pdbcode,protein_pdbcode),
                     './v2013-core/%s/%s_pocket.pdb' % (protein_pdbcode,protein_pdbcode),
                     './v2013-core/%s/%s_pocket.pdb' % (protein_pdbcode,protein_pdbcode),
                     './v2013-core/%s/%s_pocket.pdb' % (protein_pdbcode,protein_pdbcode)]
    ligand_files = ['./v2013-core/%s/%s_ligand.sdf' % (ligand_pdbcode,ligand_pdbcode),
                     './v2013-core/%s/%s_ligand.sdf' % (ligand_pdbcode,ligand_pdbcode),
                     './v2013-core/%s/%s_ligand.sdf' % (ligand_pdbcode,ligand_pdbcode),
                     './v2013-core/%s/%s_ligand.sdf' % (ligand_pdbcode,ligand_pdbcode),
                     './v2013-core/%s/%s_ligand.sdf' % (ligand_pdbcode,ligand_pdbcode),
                     './v2013-core/%s/%s_ligand.sdf' % (ligand_pdbcode,ligand_pdbcode)]    
    # load and featurize each complex
    features = self.featurizer.featurize(list(zip(ligand_files, protein_files)))
    dataset = dc.data.DiskDataset.from_numpy(features, y=labels, ids=pdbs)

    return dataset

  def _process_pdbs(self) -> Tuple[List[str], List[str], np.ndarray, List[str]]:
    data_folder = os.path.join(self.data_dir, 'v2013-core')
    index_labels_file = os.path.join(data_folder, 'pdbbind_v2013_core.csv')
      
    df = pd.read_csv(index_labels_file)
    pdbs = df.pdb_id.tolist()
    labels = np.array(df.label.tolist())
      

    if self.pocket:  # only load binding pocket
      protein_files = [
          os.path.join(data_folder, pdb, "%s_pocket.pdb" % pdb) for pdb in pdbs
      ]
    else:
      protein_files = [
          os.path.join(data_folder, pdb, "%s_protein.pdb" % pdb) for pdb in pdbs
      ]
    ligand_files = [
        os.path.join(data_folder, pdb, "%s_ligand.sdf" % pdb) for pdb in pdbs
    ]

    return (protein_files, ligand_files, labels, pdbs)


def load_input(
    featurizer: dc.feat.ComplexFeaturizer,
    splitter: Union[dc.splits.Splitter, str, None] = 'random',
    transformers: List[Union[TransformerGenerator, str]] = ['normalization'],
    reload: bool = True,
    data_dir: Optional[str] = None,
    save_dir: Optional[str] = None,
    protein_pdbcode: str ="",
    ligand_pdbcode: str="",
    pocket: bool = True,
    set_name: str = 'core',
    **kwargs
) -> Tuple[List[str], Tuple[Dataset, ...], List[dc.trans.Transformer]]:
  """Load PDBBind input.

  The PDBBind dataset includes experimental binding affinity data
  and structures for 4852 protein-ligand complexes from the "refined set"
  and 12800 complexes from the "general set" in PDBBind v2019 and 193
  complexes from the "core set" in PDBBind v2013.
  The refined set removes data with obvious problems
  in 3D structure, binding data, or other aspects and should therefore
  be a better starting point for docking/scoring studies. Details on
  the criteria used to construct the refined set can be found in [4]_.
  The general set does not include the refined set. The core set is
  a subset of the refined set that is not updated annually.

  Random splitting is recommended for this dataset.

  The input contains the columns below:

  - "ligand" - SDF of the molecular structure
  - "protein" - PDB of the protein structure
  - "CT_TOX" - Clinical trial results

  """

  loader = _PDBBindInputLoader(
      featurizer,
      splitter,
      transformers,
      PDBBIND_TASKS,
      data_dir,
      save_dir,
      pocket=pocket,
      set_name=set_name,
      **kwargs)
  return loader.load_dataset(loader.name,protein_pdbcode, ligand_pdbcode, reload)
