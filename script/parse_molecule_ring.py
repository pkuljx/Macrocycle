import pandas as pd
from rdkit import Chem

def get_max_ring_size(atomrings):
  return(max([ len(item) for item in atomrings]))

new_lines=['SMILES	ligandID	ligandname\t']
for idx, line in enumerate( open('ligand-expo.rcsb.org/Components-smiles-oe.smi').readlines()):
  if(line.count('\t')!=2):print(idx+1,'..>',line)
  else: new_lines.append(line)

with open('ligand-expo.rcsb.org/Components-smiles-oe-clean.smi','w+') as fw:
   fw.writelines(new_lines)

#SMILES  ligandID        ligandname
df = pd.read_csv('ligand-expo.rcsb.org/Components-smiles-oe-clean.smi',dtype=str,sep='\t')

exceptions=['08T', '0H2', '0JC', '0OD', '0TN', '0UE', '10R', '10S', '188', '1CU', '1FH', '1G1', '1G6', '1LN', '1MK', '1PT', '1WT', '25X', '25Y', '26E', '2AS', '2FJ', '2FK', '2GO', '2J0', '2OF', '2PT', '34B', '39B', '39E', '3JI', '3OF', '3T3', '3UQ', '3WB', '3ZZ', '4A6', '4EX', '4EY', '4HE', '4IR', '4LA', '50L', '522', '543', '5IR', '68G', '6CO', '6CQ', '72B', '76R', '7BU', '83L', '89R', '8M0', '8QJ', '8TH', '8VV', '8WV', '8ZR', '9CO', '9D7', '9JA', '9JJ', '9JM', '9OR', '9TH', '9UK', 'A71', 'A72', 'AC9', 'ALB', 'AOH', 'AUJ', 'B12', 'B1M', 'B1Z', 'B8B', 'BAZ', 'BBQ', 'BCB', 'BCL', 'BJ5', 'BOZ', 'BPT', 'BTF', 'BU0', 'BVA', 'BVR', 'BW9', 'BWU', 'C2C', 'C6J', 'C6M', 'C7P', 'CB5', 'CBY', 'CCH', 'CD1', 'CD3', 'CD5', 'CFN', 'CGO', 'CHL', 'CL0', 'CL1', 'CL2', 'CL7', 'CLA', 'CLN', 'CLZ', 'CN1', 'CNC', 'CNF', 'CO5', 'COB', 'COH', 'CON', 'COY', 'CPO', 'CPT', 'CU6', 'CUF', 'CUP', 'CUS', 'CV0', 'CWO', 'D0X', 'D6N', 'DAE', 'DAQ', 'DDH', 'DE4', 'DEF', 'DF7', 'DGQ', 'DHE', 'DJ4', 'DKE', 'DOS', 'DRU', 'DVT', 'DVW', 'DW1', 'DW2', 'DW5', 'E07', 'E52', 'E9A', 'E9D', 'EAQ', 'EJ2', 'EL9', 'ELJ', 'EQ1', 'EQQ', 'F0L', 'F0X', 'F43', 'F6C', 'FCE', 'FCI', 'FDC', 'FDE', 'FEL', 'FEM', 'FLL', 'FNE', 'FO4', 'FS2', 'FV2', 'G2O', 'G9R', 'GB0', 'GBF', 'GCR', 'GIX', 'GOV', 'GS0', 'GXW', 'GXZ', 'H1T', 'HAS', 'HB1', 'HC0', 'HC1', 'HCN', 'HE6', 'HEA', 'HEB', 'HEG', 'HEM', 'HEV', 'HF3', 'HF5', 'HIF', 'HKL', 'HNI', 'HP5', 'HUJ', 'HWS', 'I2A', 'ICE', 'ICG', 'ICH', 'ICS', 'ICZ', 'IME', 'IMF', 'IRI', 'ISW', 'J0K', 'J0N', 'J1R', 'J1S', 'J4D', 'J7T', 'J85', 'J8B', 'J8E', 'JCT', 'JGH', 'JM1', 'JQZ', 'JR3', 'JR8', 'K6G', 'KBW', 'KCO', 'KEG', 'KHK', 'KHN', 'KK5', 'KKE', 'KKH', 'KO4', 'KQB', 'KYS', 'KYT', 'L0C', 'L0E', 'L2D', 'L2M', 'L4D', 'LOS', 'LPT', 'LRU', 'LVQ', 'M10', 'M43', 'M6O', 'M7E', 'MD9', 'ME3', 'MH0', 'MI9', 'MLX', 'MM1', 'MM2', 'MM5', 'MN5', 'MN6', 'MNQ', 'MO1', 'MO2', 'MO3', 'MO4', 'MO5', 'MO6', 'MO7', 'MP1', 'MW1', 'MW2', 'MW3', 'MYW', 'N1B', 'N2N', 'N2R', 'N2S', 'N2W', 'NA2', 'NA5', 'NA6', 'NAO', 'NAW', 'NCO', 'NE5', 'NFC', 'NFV', 'NI1', 'NI2', 'NI3', 'NIK', 'NMQ', 'NRU', 'NTE', 'NUF', 'NXC', 'O1N', 'O4M', 'OC1', 'OC2', 'OC3', 'OC4', 'OC5', 'OC6', 'OC7', 'OC8', 'OCL', 'OCM', 'OCN', 'OCO', 'OEC', 'OER', 'OEX', 'OEY', 'OF1', 'OF3', 'OL3', 'OL4', 'OL5', 'OLS', 'ONP', 'OS1', 'OSW', 'OT1', 'OTQ', 'OWK', 'OXV', 'OY5', 'OY8', 'OZN', 'P5F', 'P5T', 'P6D', 'P6Q', 'P7H', 'P7Z', 'P82', 'P8B', 'P9G', 'PCD', 'PFC', 'PHF', 'PMR', 'PN8', 'PNQ', 'POR', 'PQ3', 'PQJ', 'PT7', 'PT9', 'PTE', 'PTN', 'Q2Z', 'Q38', 'Q3E', 'Q3H', 'Q3K', 'Q3N', 'Q3Q', 'Q3T', 'Q3W', 'Q41', 'Q4B', 'Q7V', 'QEB', 'QOJ', 'QQA', 'QT4', 'R1N', 'R1Z', 'RBN', 'RBU', 'RCS', 'RCZ', 'REI', 'REJ', 'REP', 'REQ', 'RHD', 'RHL', 'RHX', 'RIR', 'RKF', 'RKL', 'RKM', 'RKP', 'RLB', 'RPS', 'RR2', 'RTB', 'RTC', 'RU0', 'RU7', 'RU8', 'RUC', 'RUH', 'RUI', 'RUL', 'RUO', 'RUR', 'RUX', 'S18', 'S31', 'S32', 'S5T', 'SI0', 'SI4', 'SI7', 'SI8', 'SI9', 'SIW', 'SNF', 'SQ1', 'TBR', 'TEW', 'TIL', 'TPT', 'U8G', 'UFE', 'UTX', 'V22', 'V4A', 'V7O', 'V9G', 'VA3', 'VER', 'VFY', 'VOV', 'VPC', 'WCO', 'WJS', 'WNI', 'WO2', 'WO3', 'WUP', 'WVP', 'WXP', 'WYP', 'X33', 'X3P', 'XC3', 'XCO', 'XCU', 'XZD', 'YBT', 'YH', 'YXX', 'YXZ', 'ZAS', 'ZEM', 'ZKG', 'ZN3', 'ZNA', 'ZNB', 'ZND', 'ZNH', 'ZNO', 'ZO3', 'ZPT', 'ZRW']
ids = df['ligandID'].values

lines = ['ligandID,SMILES,maxring,numring\n']
for idx, smi in enumerate(df['SMILES'].values.astype(str)):
  if(smi == 'nan'):continue
  ligand = ids[idx]
  mol=Chem.MolFromSmiles(smi,sanitize=True)
  if mol is None and  smi != ' ' and ligand not in exceptions: 
    #exceptions.append(ligand)
    stop
  else:
    if(ligand not in exceptions): 
      rings=mol.GetRingInfo()
      atomrings = rings.AtomRings()
      if len((atomrings))==0:continue
      max_size = get_max_ring_size(atomrings)
      if(max_size>=8):
        lines.append(','.join([ligand,smi,str(max_size),str(len(atomrings))])+'\n')
        print(ligand,smi,max_size,len(atomrings))


with open('microcycle-ligands.csv','w+') as fw:
  fw.writelines(lines)
      
