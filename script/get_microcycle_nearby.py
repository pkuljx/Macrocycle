from Bio.PDB import *
from Bio.PDB.ResidueDepth import get_surface
import numpy as np
from Bio.PDB.PDBIO import Select
from Bio.PDB.ResidueDepth import min_dist
import pandas as pd
import os


class Nearbyselect:
   def __init__(self,chain_id=None):
      self.chain_id = chain_id

   def __repr__(self):
        """Represent the output as a string for debugging."""
        return "<Select all>"

   def accept_residue(self, residue):
       if(0):
        try: 
         ca = residue['CA']
         dist = min_dist(ca.get_coord(), surface)
        except:
         coords = np.asarray([a.coord for a in residue.child_list], dtype=np.float32)
         coord_center = np.average(coords, axis=0, weights=None)
         dist = min_dist(coord_center, surface)
       coords = np.asarray([a.coord for a in residue.child_list], dtype=np.float32)
       coord_center = np.average(coords, axis=0, weights=None)
       dist = min_dist(coord_center, surface)
       if(dist <=10):return True
       else: return False
     
   def accept_chain(self, chain):
       if self.chain_id is None: return 1
       else:
         if chain.get_id()==self.chain_id: return True
         else: return False


   def accept_model(self, model):
       if(model.get_id()==0): return True
       else: return False


   def accept_atom(self, atom):
        """Overload this to reject atoms for output."""
        return 1


def has_ligand(chain,ligand):
  for key in chain.child_dict.keys():
    if('H_'+ligand) in key[0]: return(key)
  return(False)


parser = PDBParser(PERMISSIVE = True, QUIET = True)

cif_parser = MMCIFParser(QUIET = True)
surface = None

io = PDBIO()
cif_io = MMCIFIO()

#ligandID,SMILES,maxring,numring,PDBs,PDBs_num
df= pd.read_csv('microcycle-ligands-pdbs.csv')
ligand_pdbs_dict=dict(zip(df['ligandID'].astype(str),df['PDBs'].astype(str)))

header_names = ['name', 'head', 'idcode', 'deposition_date',  'structure_method', 'resolution',    'compound', 'source', 'has_missing_residues',  'keywords', 'journal']

pdbl = PDBList()
result_lines = ["\t".join("ligand,pdb,chain_id,protein_name,header.head,resolution,header.name,deposition_date,structure_method,nearby_file".split(","))+"\n"]
for ligand, pdbs in ligand_pdbs_dict.items():
  pdb_list = pdbs.split()
  for pdb in pdb_list:

     #ligand = "DI0" 
     #pdb = "6xza"
     if(not os.path.exists("/data/PDBdata/pdb/pdbstructure/{}.pdb".format(pdb))):

       cif_dict=MMCIF2Dict.MMCIF2Dict("./pdb/{}.cif".format(pdb))
       serial_entity_dict = dict(zip(cif_dict['_atom_site.id'], cif_dict['_atom_site.label_entity_id']))
       
       #pdbl.retrieve_pdb_file(pdb,pdir='./pdb', file_format='pdb')

#_entity.id >>> ['1', '2', '3', '4', '5']
#_entity.type >>> ['polymer', 'polymer', 'non-polymer', 'non-polymer', 'non-polymer']
#_entity.src_method >>> ['man', 'man', 'syn', 'syn', 'syn']
#_entity.pdbx_description >>> ['Myosin-A', 'Actin-1', "ADENOSINE-5'-DIPHOSPHATE", 'MAGNESIUM ION', 'Jasplakinolide']
       entity_protein_dict = dict(zip([str(item) for item in cif_dict['_entity.id']],cif_dict['_entity.pdbx_description']))
       print(entity_protein_dict)
       structure = cif_parser.get_structure(pdb,"./pdb/{}.cif".format(pdb))

       model = structure[0]
       header = structure.header
       def get_chain_protein(chain):

          for atom in chain.get_atoms():
             print([atom.get_id()])
             if(atom.get_id() in ["CA","P"]):
              print(atom)
              print(chain.get_id(),atom.get_serial_number(),str(serial_entity_dict[str(atom.get_serial_number())]))
              return(entity_protein_dict[str(serial_entity_dict[str(atom.get_serial_number())])])
  
  
       
       for chain in model.get_chains():

            if chain.get_id() !="A2":continue
            #print(chain.child_list)
            #print(chain.child_dict)
            #print(chain.get_residues)
            protein_name =get_chain_protein(chain)
            chain_id = chain.get_id()
            ligand_id = has_ligand(chain, ligand)
            if(not ligand_id):continue
            print(chain_id)
            print([ligand,pdb,chain_id,protein_name,str(header['head']),str(header["resolution"]),header["name"],str(header["deposition_date"]),header["structure_method"],"{}_{}_{}.cif".format(ligand,pdb,chain_id)])
            print(chain_id)
            result_lines.append('\t'.join([ligand,pdb,chain_id,protein_name,str(header['head']),str(header["resolution"]),header["name"],str(header["deposition_date"]),header["structure_method"],"{}_{}_{}.cif".format(ligand,pdb,chain_id)])+"\n")
            print(chain_id)
            print('\t'.join([ligand,pdb,chain_id,protein_name,str(header['head']),str(header["resolution"]),header["name"],str(header["deposition_date"]),header["structure_method"],"{}_{}_{}.cif".format(ligand,pdb,chain_id)]))
            residue = chain[ligand_id]
            surface = np.asarray([a.coord for a in residue.child_list], dtype=np.float32)
            cif_io.set_structure(structure)
            cif_io.save("./ligand_nearby/{}_{}_{}.cif".format(ligand,pdb,chain_id), Nearbyselect(chain_id=None))
       #break     
     else:    
       
       structure = parser.get_structure(pdb,"/data/PDBdata/pdb/pdbstructure/{}.pdb".format(pdb))
       #ligand = '5I1'
       #pdb = '5acb'
       #structure = parser.get_structure(pdb,"5ACB.pdb") 
       model = structure[0]
       header = structure.header
 
       chain_protain_dict={}
       for value in header['compound'].values():
         for c in value["chain"].upper().split(","):
          chain_protain_dict[c.strip()] = value['molecule']
    

       for chain in model.get_chains():
            #print(chain.child_list)
            #print(chain.child_dict)
            #print(chain.get_residues)

            chain_id = chain.get_id()
            #1fq5
            if chain_id not in chain_protain_dict:continue
            protein_name = chain_protain_dict[chain_id]
            ligand_id = has_ligand(chain, ligand)
            if(not ligand_id):continue
            result_lines.append('\t'.join([ligand,pdb,chain_id,protein_name,header['head'],str(header["resolution"]),header["name"],str(header["deposition_date"]),header["structure_method"],"{}_{}_{}.pdb".format(ligand,pdb,chain_id)])+"\n")
           
            print('\t'.join([ligand,pdb,chain_id,protein_name,header['head'],str(header["resolution"]),header["name"],str(header["deposition_date"]),header["structure_method"],"{}_{}_{}.pdb".format(ligand,pdb,chain_id)]))
            residue = chain[ligand_id]
            surface = np.asarray([a.coord for a in residue.child_list], dtype=np.float32)
            io.set_structure(structure)
            io.save("./ligand_nearby/{}_{}_{}.pdb".format(ligand,pdb,chain_id), Nearbyselect(chain_id=None))    
       #break  



       


with open("microcycle-ligands-pdbs-nearby.tsv",'w+') as fw:
  fw.writelines(result_lines)

