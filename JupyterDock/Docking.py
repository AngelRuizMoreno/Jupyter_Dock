from pymol import cmd
from openbabel import pybel

from rdkit import Chem

import os, shutil, subprocess

from .DockingHelpers import dok_to_sdf

from importlib_resources import files
SMINA = files("JupyterDock.bin").joinpath("smina")
LEDOCK = files("JupyterDock.bin").joinpath("ledock_linux_x86")

class Docker():
    def __init__(self,receptor_file:str,ligands_file:str,selection:str,extending:float=5.0):
        self.receptor_file =receptor_file
        self.ligands_file = ligands_file
        self.selection = selection
        self.extending = extending
        
        self.docking_box={}

        cmd.load(filename=self.receptor_file,format='pdb',object='receptor')
        cmd.load(filename=self.ligands_file,format='mol2',object='lig')

        ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent(self.selection)

        minX = minX - float(self.extending)
        minY = minY - float(self.extending)
        minZ = minZ - float(self.extending)
        maxX = maxX + float(self.extending)
        maxY = maxY + float(self.extending)
        maxZ = maxZ + float(self.extending)

        SizeX = maxX - minX
        SizeY = maxY - minY
        SizeZ = maxZ - minZ
        CenterX =  (maxX + minX)/2
        CenterY =  (maxY + minY)/2
        CenterZ =  (maxZ + minZ)/2

        cmd.delete('all')

        self.docking_box['smina']={'--center_x':CenterX,'--center_y': CenterY, '--center_z': CenterZ, '--size_x':SizeX,'--size_y': SizeY,'--size_z': SizeZ}
        self.docking_box['ledock']= {'minX':minX, 'maxX': maxX,'minY':minY, 'maxY':maxY,'minZ':minZ,'maxZ':maxZ}

        
    def run_smina(self,params:dict={'--exhaustiveness':8,'--num_modes': 20},verbose:bool=False):
        
        self.smina_params={'-r' : self.receptor_file,
                           '-l' : self.ligands_file,
                           '-o': f"{self.ligands_file.split('.')[0]}_smina_results.sdf"}
        
        self.smina_params.update(self.docking_box['smina'])
        
        self.smina_params.update(params)
            
        cmd = [] 
        for k, v in self.smina_params.items():
            if k != '' or v!='':
                cmd.append(str(k)) 
                cmd.append(str(v))

        self.result = subprocess.run([SMINA]+cmd, capture_output=True-verbose, text=True)
                
    def run_ledock(self,rmsd:float=1.0, num_modes:int=20):
        self.ledock_params={'Receptor' : self.receptor_file,
                            'Ligands' : self.ligands_file,
                            'Results': f"{self.ligands_file.split('.')[0]}_ledock_results.sdf",
                            'RMSD' : str(rmsd),
                            'num_modes' : str(num_modes),
                           }
        
        self.ledock_params.update(self.docking_box['ledock'])
        
        
        
        tmp_dir=os.path.join('.','_dock_tmp')
        if os.path.isdir(tmp_dir) ==False:
            os.mkdir(tmp_dir)
        else:
            pass

        
        mols = [m for m in pybel.readfile(filename=self.ligands_file,format=self.ligands_file.split('.')[1])]
        
        for i,mol in enumerate(mols):
            if mol.title:
                out=pybel.Outputfile(filename=os.path.join(tmp_dir,f"{mol.title}.mol2"),format='mol2',overwrite=True)
            else:
                out=pybel.Outputfile(filename=os.path.join(tmp_dir,f"ligand_{i+1}.mol2"),format='mol2',overwrite=True)
            out.write(mol)
            out.close()
        
        
        with open(os.path.join(tmp_dir,'ligands.lst'),'w') as lig_lst:
            for file in os.listdir(tmp_dir):
                if 'mol2' in file:
                    lig_lst.write(os.path.join(tmp_dir,file))
        lig_lst.close()
        
        
        file=f'''Receptor
{self.ledock_params['Receptor']}

RMSD
{self.ledock_params['RMSD']}

Binding pocket
{self.ledock_params['minX']} {self.ledock_params['maxX']}
{self.ledock_params['minY']} {self.ledock_params['maxY']}
{self.ledock_params['minZ']} {self.ledock_params['maxZ']}

Number of binding poses
{self.ledock_params['num_modes']}

Ligands list
{os.path.join(tmp_dir,'ligands.lst')}

END'''
        
        with open(os.path.join(tmp_dir,'dock.in'),'w') as out:
            out.write(file)
        
        out.close()
            
        
        result = subprocess.run([LEDOCK, os.path.join(tmp_dir,'dock.in')], capture_output=True, text=True)
        
        
        for dok_file in os.listdir(tmp_dir):
            if 'dok' in dok_file:
                dok_to_sdf(os.path.join(tmp_dir,dok_file),f"{self.ligands_file.split('.')[0]}_ledock_{dok_file.replace('dok','sdf')}")
        
        shutil.rmtree(tmp_dir)