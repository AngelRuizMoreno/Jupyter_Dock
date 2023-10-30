import py3Dmol
from rdkit import Chem


class VisualizeDocking():
    
    def __init__(self,receptor:str, ligands:str,reference:str=None, pose:int = 1):
        
        self.receptor = receptor
        self.ligands = ligands
        self.reference = reference
        self.pose = pose
        
    def view (self):

        viewer = py3Dmol.view()
        viewer.removeAllModels()
        viewer.setViewStyle({'style':'outline','color':'black','width':0.1})

        viewer.addModel(open(self.receptor,'r').read(),format='pdb')
        
        protein=viewer.getModel()
        protein.setStyle({'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}})
        
        viewer.addSurface(py3Dmol.VDW,{'opacity':0.6,'color':'white'})
        
        if isinstance(self.reference,str):
    
            viewer.addModel(open(self.reference,'r').read(),format='mol2')

            mol_ref = viewer.getModel()
            mol_ref.setStyle({},{'stick':{'colorscheme':'redCarbon','radius':0.2}})
        
        else:
            pass

        
        mol_dock=Chem.SDMolSupplier(self.ligands)

        pose=Chem.MolToMolBlock(mol_dock[self.pose-1],False)

        

        viewer.addModel(pose,'mol')
        
        p = viewer.getModel()
        p.setStyle({},{'stick':{'colorscheme':'cyanCarbon','radius':0.2}})
        
        print(f"Docking_result : Cyan | Pose Number: {mol_dock[self.pose-1].GetProp('Pose')} | Score : {mol_dock[self.pose-1].GetProp('Score')}")
        
        viewer.zoomTo()
        viewer.show()
        
        return viewer