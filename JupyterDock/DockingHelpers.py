import py3Dmol
from pymol import cmd
from openbabel import pybel

from rdkit import Chem
from rdkit.Chem import AllChem,rdFMCS, rdMolAlign

from pdbfixer import PDBFixer
from openmm.app import PDBFile

import MDAnalysis as mda
from MDAnalysis.coordinates import PDB

import random, math

import numpy as np

__all__=['get_docking_box','sanitize_protein','generate_ledock_file','dok_to_sdf', 'pdbqt_to_sdf','compute_inplace_rmsd', 'get_scaffold_based_conformers']


def dok_to_sdf (dok_file:str,output:str):
    '''
    Description: The function dok_to_sdf converts a DOK file to an SDF file.

    Parameters:

    - dok_file: A string that specifies the name of the input DOK file. A DOK file is a text file that contains the coordinates and properties of a set of molecules in PDB format. If None, the function will raise an exception.
    - output: A string that specifies the name of the output SDF file. An SDF file is a text file that contains the structure and data of a set of molecules in a standard format. If None, the function will use the same name as the input file with the ‘.sdf’ extension.
    - Return: The function does not return anything, but it writes the output SDF file to the disk.

    Example: To convert a DOK file named ‘example.dok’ to an SDF file named ‘example.sdf’, you can call the function as follows:

    >>> dok_to_sdf(dok_file='example.dok', output='example.sdf')
    '''

    out=pybel.Outputfile(filename=output,format='sdf',overwrite=True)

    with open(dok_file, 'r') as f:
        doc=f.readlines()

    start=[index for (index,p) in enumerate(doc) if 'REMARK Cluster' in p]
    end=[index-1 for (index,p) in enumerate(doc) if 'END' in p]
    

    interval=list(zip(start,end))
    for pair in interval:
        block = ",".join(doc[pair[0]:pair[1]]).replace(',','')
        m=pybel.readstring(format='pdb',string=block)
        m.title=dok_file.split('/')[0].split('.')[0]
        m.data.update({'minimizedAffinity':m.data['REMARK'].split()[6]})
        del m.data['REMARK']

        out.write(m)

    out.close()
  

def pdbqt_to_sdf(pdbqt_file:str,output:str):
    
    '''
    The function pdbqt_to_sdf converts a PDBQT file to an SDF file.

    Parameters:
    - pdbqt_file: A string that specifies the name of the input PDBQT file. A PDBQT file is a text file that contains the coordinates and properties of a set of molecules in a format suitable for AutoDock Vina. If None, the function will raise an exception.
    - output: A string that specifies the name of the output SDF file. An SDF file is a text file that contains the structure and data of a set of molecules in a standard format. If None, the function will use the same name as the input file with the ‘.sdf’ extension.
    
    Return: The function does not return anything, but it writes the output SDF file to the disk.

    Example: To convert a PDBQT file named ‘example.pdbqt’ to an SDF file named ‘example.sdf’, you can call the function as follows:

    >>> pdbqt_to_sdf(pdbqt_file='example.pdbqt', output='example.sdf')
    '''

    results = [m for m in pybel.readfile(filename=pdbqt_file,format='pdbqt')]
    out=pybel.Outputfile(filename=output,format='sdf',overwrite=True)
    for pose in results:

        pose.data.update({'Pose':pose.data['MODEL']})
        pose.data.update({'Score':pose.data['REMARK'].split()[2]})
        del pose.data['MODEL'], pose.data['REMARK'], pose.data['TORSDO']

        out.write(pose)
    out.close()


def compute_inplace_rmsd (ref:Chem.Mol,target:Chem.Mol):
    
    r=rdFMCS.FindMCS([ref,target])
    
    a=ref.GetSubstructMatch(Chem.MolFromSmarts(r.smartsString))
    b=target.GetSubstructMatch(Chem.MolFromSmarts(r.smartsString))   
    amap=list(zip(a,b))
    
    distances=[]
    for atomA, atomB in amap:
        pos_A=ref.GetConformer().GetAtomPosition (atomA)
        pos_B=target.GetConformer().GetAtomPosition (atomB)
        coord_A=np.array((pos_A.x,pos_A.y,pos_A.z))
        coord_B=np.array ((pos_B.x,pos_B.y,pos_B.z))
        dist_numpy = np.linalg.norm(coord_A-coord_B)        
        distances.append(dist_numpy)
        
    rmsd=math.sqrt(1/len(distances)*sum([i*i for i in distances]))
    
    return rmsd

def get_scaffold_based_conformers(smiles:str, anchor:Chem.Mol, num_confs:int, output:str, rmsd_threshold:int=0.75):
    mol = Chem.MolFromSmiles(smiles,sanitize=True)
    mol=Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)

    constrain = Chem.SDMolSupplier(anchor,sanitize=True)[0]

    r = rdFMCS.FindMCS([mol, constrain])
    a = mol.GetSubstructMatch(Chem.MolFromSmarts(r.smartsString))
    b = constrain.GetSubstructMatch(Chem.MolFromSmarts(r.smartsString))
    amap = list(zip(a, b))
    coors = dict()

    for a,b in amap:
        coors[a] = constrain.GetConformer().GetAtomPosition(b)

    w = Chem.SDWriter(output)

    mol.UpdatePropertyCache()
    constrain.UpdatePropertyCache()

    confs = AllChem.EmbedMultipleConfs(mol,
        numConfs=int(num_confs),
        coordMap=coors,
        pruneRmsThresh=0.75,
        useExpTorsionAnglePrefs=True,
        useBasicKnowledge=True)

    for element in confs:
        Chem.SanitizeMol(mol)
        rmsd = AllChem.GetBestRMS(mol,constrain,element,0,map=[list(amap)])
        if rmsd<=float(rmsd_threshold):
            w.write(mol, confId=element)
    w.close()

'''
def get_3D_view (receptor_file='',rec_opts={'format':'pdb'},docking_results='',refMol='',refMol_opts={'format':'mol2'},pose=[0]):

    view = py3Dmol.view()
    view.removeAllModels()
    view.setViewStyle({'style':'outline','color':'black','width':0.1})

    view.addModel(open(receptor_file,'r').read(),**rec_opts)
    Prot=view.getModel()
    Prot.setStyle({'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}})
    view.addSurface(py3Dmol.VDW,{'opacity':0.6,'color':'white'})

    if refMol:
        view.addModel(open(refMol,'r').read(),**refMol_opts)
        ref_m = view.getModel()
        ref_m.setStyle({},{'stick':{'colorscheme':'greenCarbon','radius':0.2}})

    if pose:
        results=Chem.SDMolSupplier(docking_results)
        for index in pose:

            color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])]    
            p=Chem.MolToMolBlock(results[index])
            #print (results[index].GetProp('REMARK'))

            view.addModel(p,'mol')
            x = view.getModel()
            x.setStyle({},{'stick':{'color':color[0],'radius':0.1}})

    view.zoomTo()
    view.show()
'''