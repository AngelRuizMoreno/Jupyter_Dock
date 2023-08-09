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

def get_docking_box(selection='sele', extending = 6.0, software='vina'):
    """
    Calculates the docking box for a given pymol selection: atom(s), residue(s), molecule(s).

    Parameters
    ----------
    selection : str, optional
        The name of the selection of molecules in PyMOL. The default is 'sele'.
    extending : float, optional
        The amount to extend the docking box from the minimum and maximum coordinates of the selection. The default is 6.0.
    software : str, optional
        The software to use for docking. The options are 'vina', 'ledock', or 'both'. The default is 'vina'.

    Returns
    -------
    dict or tuple
        The format of the docking box depending on the software parameter.
        If software is 'vina', it returns a dictionary with keys 'center_x', 'center_y', 'center_z' for the center coordinates, and another dictionary with keys 'size_x', 'size_y', 'size_z' for the size of the docking box.
        If software is 'ledock', it returns three tuples, each containing two values for the minimum and maximum coordinates of the docking box in x, y, and z axes, respectively.
        If software is 'both', it returns both formats as a tuple of two elements.
        If software is not one of the options, it prints an error message and returns None.

    Example
    -------
    >>> get_docking_box(selection='ligand', extending=4.0, software='both')
    (({'center_x': 10.5, 'center_y': 20.0, 'center_z': 30.5}, {'size_x': 12.0, 'size_y': 16.0, 'size_z': 20.0}), ({'minX': 4.5, 'maxX': 16.5}, {'minY': 12.0, 'maxY': 28.0}, {'minZ': 20.5, 'maxZ': 40.5}))
    """
    
    ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent(selection)

    minX = minX - float(extending)
    minY = minY - float(extending)
    minZ = minZ - float(extending)
    maxX = maxX + float(extending)
    maxY = maxY + float(extending)
    maxZ = maxZ + float(extending)
    
    SizeX = maxX - minX
    SizeY = maxY - minY
    SizeZ = maxZ - minZ
    CenterX =  (maxX + minX)/2
    CenterY =  (maxY + minY)/2
    CenterZ =  (maxZ + minZ)/2
    
    cmd.delete('all')
    
    if software == 'vina':
        return {'center_x':CenterX,'center_y': CenterY, 'center_z': CenterZ},{'size_x':SizeX,'size_y': SizeY,'size_z': SizeZ}
    elif software == 'ledock':
        return {'minX':minX, 'maxX': maxX},{'minY':minY, 'maxY':maxY}, {'minZ':minZ,'maxZ':maxZ}
    elif software == 'both':
        return ({'center_x':CenterX,'center_y': CenterY, 'center_z': CenterZ},{'size_x':SizeX,'size_y': SizeY,'size_z': SizeZ}),({'minX':minX, 'maxX': maxX},{'minY':minY, 'maxY':maxY}, {'minZ':minZ,'maxZ':maxZ})
    
    else:
        print('software options must be "vina", "ledock" or "both"')


def sanitize_protein(filename='',addHs_pH=7.4,output='',try_renumberResidues=False):
    """
    Sanitizes a protein structure file by using PDBFixer and MDAnalysis.

    Parameters
    ----------
    filename : str
        The name of the input file in PDB format.
    addHs_pH : float, optional
        The pH value for adding missing hydrogens to the protein. The default is 7.4.
    output : str
        The name of the output file in PDB format.
    try_renumberResidues : bool, optional
        Whether to try to renumber the residues in the output file according to the original file. The default is False.

    Returns
    -------
    None

    Example
    -------
    >>> sanitize_protein(filename='1a2b.pdb', output='1a2b_fixed.pdb', try_renumberResidues=True)
    """

    fix = PDBFixer(filename=filename)
    fix.findMissingResidues()
    fix.findNonstandardResidues()
    fix.replaceNonstandardResidues()
    fix.removeHeterogens(True)
    fix.findMissingAtoms()
    fix.addMissingAtoms()
    fix.addMissingHydrogens(addHs_pH)
    PDBFile.writeFile(fix.topology, fix.positions, open(output, 'w'))

    if try_renumberResidues == True:
        try:
            original=mda.Universe(filename)
            from_fix=mda.Universe(output)

            resNum=[res.resid for res in original.residues]
            for idx,res in enumerate(from_fix.residues):
                res.resid = resNum[idx]

            save=PDB.PDBWriter(filename=output)
            save.write(from_fix)
            save.close()
        except Exception:
            print('Not possible to renumber residues, check exception for extra details')
        

def generate_ledock_file(receptor:str='pro.pdb',rmsd:float=1.0,x:list=[0,0],y:list=[0,0],z:list=[0,0], n_poses:int=10, l_list:list=[],l_list_outfile:str='',out:str='dock.in'):
     """
    Generates an input file for LeDock docking program.

    Parameters
    ----------
    receptor : str
        The name of the receptor file in PDB format. The default is 'pro.pdb'.
    rmsd : float, optional
        The root mean square deviation threshold for clustering the docking poses. The default is 1.0.
    x : list, optional
        The minimum and maximum coordinates of the binding pocket in x axis. The default is [0,0].
    y : list, optional
        The minimum and maximum coordinates of the binding pocket in y axis. The default is [0,0].
    z : list, optional
        The minimum and maximum coordinates of the binding pocket in z axis. The default is [0,0].
    n_poses : int, optional
        The number of binding poses to output. The default is 10.
    l_list : list
        The names of the ligand files in MOL2 format.
    l_list_outfile : str
        The name of the output file that contains the ligand list.
    out : str, optional
        The name of the output file in LEDock format. The default is 'dock.in'.

    Returns
    -------
    None

    Example
    -------
    >>> generate_ledock_file(receptor='protein.pdb', rmsd=2.0, x=[10,20], y=[15,25], z=[5,15], n_poses=5, l_list=['ligand1.mol2', 'ligand2.mol2'], l_list_outfile='ligands.txt', out='ledock.in')
    """
        
    rmsd=str(rmsd)
    x=[str(x) for x in x]
    y=[str(y) for y in y]
    z=[str(z) for z in z]
    n_poses=str(n_poses)

    with open(l_list_outfile,'w') as l_out:
        for element in l_list:
            l_out.write(element)
    l_out.close()

    file=[
        'Receptor\n',
        receptor + '\n\n',
        'RMSD\n',
        rmsd +'\n\n',
        'Binding pocket\n',
        x[0],' ',x[1],'\n',
        y[0],' ',y[1],'\n',
        z[0],' ',z[1],'\n\n',
        'Number of binding poses\n',
        n_poses + '\n\n',
        'Ligands list\n',
        l_list_outfile + '\n\n',
        'END']
    
    with open(out,'w') as output:
        for line in file:
            output.write(line)
    output.close()


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
        doc=[line for line in f.readlines()]
    
    doc=[line.replace(line.split()[2],line.split()[2].upper()) if 'ATOM' in line else line for line in doc]
    
    start=[index for (index,p) in enumerate(doc) if 'REMARK Cluster' in p]
    finish=[index-1 for (index,p) in enumerate(doc) if 'REMARK Cluster' in p]
    finish.append(len(doc))

    interval=list(zip(start,finish[1:]))
    for num,i in enumerate(interval):
        block = ",".join(doc[i[0]:i[1]]).replace(',','')

        m=pybel.readstring(format='pdb',string=block)
        
        m.data.update({'Pose':m.data['REMARK'].split()[4]})
        m.data.update({'Score':m.data['REMARK'].split()[6]})
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