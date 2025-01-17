import pandas as pd



ligands_not_in_pdbbank={'PI2', 'FRR', 'PI3', 'TXG', '3L2', 'QB4', '48D', 'M2D', '4CQ', '0EZ', '08L', 'BH1', 'O8K', '08X', 'CNG', 'RPX', 'PI1', '53P', '3H3', 'GMY', 'AQJ', 'FE8', 'FWT', 'ECT', 'L2A'}

df_ligands = pd.read_csv('microcycle-ligands.csv')
df_ligands=df_ligands[~df_ligands['ligandID'].isin(ligands_not_in_pdbbank)]
df_ligands.set_index('ligandID',inplace=True)
ligands  = list(df_ligands.index)
print(ligands)
print(df_ligands.loc['ZPN'].to_dict())
#ligandID        PDBs

ligand_pdbs_dict = {}
df_ligand_pdbs = pd.read_csv('ligand-expo.rcsb.org/cc-to-pdb.tdd',sep='\t')
for ligand, pdbs in zip(df_ligand_pdbs['ligandID'],df_ligand_pdbs['PDBs']):
  if(ligand not in ligands): continue
  print(ligand,pdbs.split())
  ligand_pdbs_dict[ligand]=pdbs.split()

print(set(ligands)-set(ligand_pdbs_dict.keys()))
ligand_pdbs = []
ligand_pdbs_num =[]
for ligand in ligands:
  ligand_pdbs.append(' '.join(ligand_pdbs_dict[ligand]))
  ligand_pdbs_num.append(len(ligand_pdbs_dict[ligand]))
df_ligands['PDBs'] = ligand_pdbs

df_ligands['PDBs_num'] = ligand_pdbs_num
df_ligands.to_csv('microcycle-ligands-pdbs.csv')
