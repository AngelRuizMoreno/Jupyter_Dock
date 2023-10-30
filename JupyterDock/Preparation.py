from pymol import cmd

from pdbfixer import PDBFixer
from openmm.app import PDBFile
import MDAnalysis as mda
from MDAnalysis.coordinates import PDB

from openbabel import pybel

from importlib_resources import files
LEPRO = files("JupyterDock.bin").joinpath("lepro_linux_x86")



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
    
class SanitizeProtein():
    def __init__(self,protein_file:str=None):
        self.protein_file=protein_file
        self.output_file=self.protein_file.replace('.pdb','_sanitized.pdb')
    
    def lepro_sanitization(self):
        
        result = subprocess.run([LEPRO, self.protein_file,'-p'], capture_output=True, text=True)
        os.rename('pro.pdb',self.output_file)
        os.remove('dock.in')
    
    def pdbfixer_sanitization(self,pH:float=7.4):
        
        fix = PDBFixer(filename=self.protein_file)
        fix.findMissingResidues()
        fix.findNonstandardResidues()
        fix.replaceNonstandardResidues()
        fix.removeHeterogens(True)
        fix.findMissingAtoms()
        fix.addMissingAtoms()
        fix.addMissingHydrogens(pH)
        PDBFile.writeFile(fix.topology, fix.positions, open(self.output_file, 'w'))
        try:
            original=mda.Universe(self.protein_file)
            from_fix=mda.Universe(self.output_file)

            resNum=[res.resid for res in original.residues]
            for idx,res in enumerate(from_fix.residues):
                res.resid = resNum[idx]

            save=PDB.PDBWriter(filename=self.output_file)
            save.write(from_fix)
            save.close()
        except Exception:
            print('Not possible to renumber residues, check exception for extra details')
            pass
        


class SanitizeLigands():
    def __init__(self, ligands_file:str=None):
        self.ligands_file=ligands_file
        self.output_file=f"{self.ligands_file.split('.')[0]}_sanitized.mol2"
    
    def run_sanitization(self):
        
        mols= [m for m in pybel.readfile(filename=self.ligands_file,format=self.ligands_file.split('.')[1])]
        
        out=pybel.Outputfile(filename=self.output_file,format="mol2",overwrite=True)
        
        self.ligands=[]
        for i,m in enumerate(mols):
            m.convertdbonds()
            m.addh()
            if m.title == '':
                m.title = f"{self.ligands_file.split('/')[-1].split('.')[0]}_Lig{i+1}"
            else: pass
            
            out.write(m)
            self.ligands.append(m)
        out.close()