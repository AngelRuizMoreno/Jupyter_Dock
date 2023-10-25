from pymol import cmd

def fetch_protein_from_pdb(entry:str='1AZ8',extract_ligand:bool=True, assembly:int=1):
    if assembly == 0:
        cmd.fetch(code=entry)
    else:
        cmd.fetch(code=entry,type=f"pdb{assembly}")
        
    cmd.select(name='Prot',selection='polymer.protein')
    if extract_ligand == True:
        cmd.select(name='Lig',selection='organic')
        cmd.save(filename=f'{entry}_lig.mol2',format='mol2',selection='Lig')
    
    cmd.save(filename=f'{entry}_clean.pdb',format='pdb',selection='Prot')
    cmd.delete('all')