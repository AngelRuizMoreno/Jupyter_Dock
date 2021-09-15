import py3Dmol
from pymol import cmd
from openbabel import pybel

from rdkit import Chem
from rdkit.Chem import AllChem,rdFMCS, Draw

from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

import MDAnalysis as mda
from MDAnalysis.coordinates import PDB

import random

def getbox(selection='sele', extending = 6.0, software='vina'):
    
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


def fix_protein(filename='',addHs_pH=7.4,output='',try_renumberResidues=False):

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
            original=mda.Universe(output)
            from_fix=mda.Universe(output)

            resNum=[res.resid for res in original.residues]
            for idx,res in enumerate(from_fix.residues):
                res.resid = resNum[idx]

            save=PDB.PDBWriter(filename=output)
            save.write(from_fix)
            save.close()
        except Exception:
            print('SNot possible to renumber residues, check excepton for extra details')
        

def generate_ledock_file(receptor='pro.pdb',rmsd=1.0,x=[0,0],y=[0,0],z=[0,0], n_poses=10, l_list=[],l_list_outfile='',out='dock.in'):
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



def DOKMolSupplier (file=None):
    mols=[]
    with open(file, 'r') as f:
        doc=[line for line in f.readlines()]
    doc=[line.replace(line.split()[2],line.split()[2].upper()) if 'ATOM' in line else line for line in doc]
    
    scores=[line.split()[-2] for line in doc if 'REMARK Cluster' in line]
    poses=[line.split()[2] for line in doc if 'REMARK Cluster' in line]

    start=[index for (index,p) in enumerate(doc) if 'REMARK Cluster' in p]
    finish=[index-1 for (index,p) in enumerate(doc) if 'REMARK Cluster' in p]
    finish.append(len(doc))

    interval=list(zip(start,finish[1:]))
    for num,i in enumerate(interval):
        try:
            block = ",".join(doc[i[0]:i[1]]).replace(',','')

            m=pybel.readstring(format='pdb',string=block)
            mols.append(m)
        except Exception:
            pass
    return(mols)

def dok_to_sdf (dok_file=None,output=None):
    mols=DOKMolSupplier(dok_file)
    out=pybel.Outputfile(filename=output,format='sdf',overwrite=True)
    for pose in mols:
        pose.addh()
        out.write(pose)
    out.close()

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