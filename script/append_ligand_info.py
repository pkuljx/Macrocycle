from Bio.PDB import *
from Bio.PDB.ResidueDepth import get_surface
import numpy as np
from Bio.PDB.PDBIO import Select
from Bio.PDB.ResidueDepth import min_dist
import pandas as pd
import os



parser = PDBParser(PERMISSIVE = True, QUIET = True)
cif_parser = MMCIFParser(QUIET = True)

io = PDBIO()
cif_io = MMCIFIO()



#ligandID,SMILES,maxring,numring,PDBs,PDBs_num
df= pd.read_csv('microcycle-ligands-pdbs.csv')
ligand_smi_dict=dict(zip(df['ligandID'].astype(str),df['SMILES'].astype(str)))



df = pd.read_csv("microcycle-ligands-pdbs-nearby.tsv",sep='\t')

smi_list=[]
nearby_chain_num_list =[]
for nearby_file in df['nearby_file']:
  print(nearby_file)
  ligand = nearby_file.split('_')[0]
  pdb = nearby_file.split('_')[1]

  smi = ligand_smi_dict[ligand]
  if nearby_file.endswith('pdb'):
     structure = parser.get_structure(pdb,'ligand_nearby/'+nearby_file)
  if nearby_file.endswith('cif'):
     structure = cif_parser.get_structure(pdb,'ligand_nearby/'+nearby_file)
  chains=[]
  for model in structure.get_models():
    for chain in model.get_chains():
      chains.append(chain.get_id())
  print(chains)
  smi_list.append(smi)
  nearby_chain_num_list.append(len(chains))

df['SMILES'] = smi_list
df['nearby_chain_num'] = nearby_chain_num_list

df.to_csv("microcycle-ligands-pdbs-nearby_info.csv",sep=',',index=False)
