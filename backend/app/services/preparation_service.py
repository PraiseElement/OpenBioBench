"""
Preparation Service for protein and ligand preparation.
Handles cleaning, protonation, and optimization.
"""
from typing import Dict, List, Optional
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from io import StringIO
import logging

logger = logging.getLogger(__name__)


class PreparationService:
    """Service for preparing proteins and ligands for docking."""
    
    @staticmethod
    def clean_protein(pdb_content: str, options: dict = None) -> Dict:
        """
        Clean a protein structure for docking.
        
        Options:
            remove_water: Remove water molecules (default: True)
            remove_heteroatoms: Remove non-protein atoms (default: False)
            keep_ligands: List of ligand names to keep
            select_chains: List of chains to keep (None = all)
        """
        options = options or {}
        remove_water = options.get('remove_water', True)
        remove_heteroatoms = options.get('remove_heteroatoms', False)
        keep_ligands = options.get('keep_ligands', [])
        select_chains = options.get('select_chains', None)
        
        lines = pdb_content.strip().split('\n')
        cleaned_lines = []
        
        removed_waters = 0
        removed_heteroatoms = 0
        kept_atoms = 0
        
        for line in lines:
            # Keep headers and other records
            if not line.startswith('ATOM') and not line.startswith('HETATM'):
                if not (line.startswith('ANISOU')):  # Skip ANISOU records
                    cleaned_lines.append(line)
                continue
            
            # Parse chain if specified
            if len(line) > 21:
                chain = line[21]
                if select_chains and chain not in select_chains:
                    continue
            
            # Handle HETATM records
            if line.startswith('HETATM'):
                residue_name = line[17:20].strip() if len(line) > 20 else ""
                
                # Remove waters
                if remove_water and residue_name in ['HOH', 'WAT', 'H2O', 'DOD']:
                    removed_waters += 1
                    continue
                
                # Remove other heteroatoms unless in keep list
                if remove_heteroatoms and residue_name not in keep_ligands:
                    removed_heteroatoms += 1
                    continue
            
            cleaned_lines.append(line)
            kept_atoms += 1
        
        # Add END record if not present
        if cleaned_lines and not cleaned_lines[-1].startswith('END'):
            cleaned_lines.append('END')
        
        return {
            "pdb_content": '\n'.join(cleaned_lines),
            "stats": {
                "removed_waters": removed_waters,
                "removed_heteroatoms": removed_heteroatoms,
                "kept_atoms": kept_atoms
            }
        }
    
    @staticmethod
    def add_hydrogens_to_protein(pdb_content: str, ph: float = 7.0) -> Dict:
        """
        Add hydrogen atoms to protein (placeholder - uses simple approach).
        
        Production would use PDB2PQR or reduce.
        """
        # For now, just return the structure with a note
        # Real implementation would call PDB2PQR service
        lines = pdb_content.strip().split('\n')
        
        # Count existing hydrogens
        h_count = sum(1 for line in lines 
                     if (line.startswith('ATOM') or line.startswith('HETATM')) 
                     and len(line) > 77 and line[76:78].strip() == 'H')
        
        return {
            "pdb_content": pdb_content,
            "hydrogens_present": h_count,
            "ph": ph,
            "note": "Hydrogen addition requires PDB2PQR. Structure returned as-is."
        }
    
    @staticmethod
    def prepare_ligand_from_smiles(smiles: str, options: dict = None) -> Dict:
        """
        Prepare a ligand from SMILES for docking.
        
        Options:
            optimize: Run energy minimization (default: True)
            ph: Protonation pH (default: 7.0)
            num_conformers: Number of conformers to generate (default: 1)
        """
        options = options or {}
        optimize = options.get('optimize', True)
        num_conformers = options.get('num_conformers', 1)
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates
        try:
            if num_conformers > 1:
                AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, randomSeed=42)
            else:
                AllChem.EmbedMolecule(mol, randomSeed=42)
        except Exception as e:
            raise ValueError(f"Failed to generate 3D coordinates: {e}")
        
        # Optimize geometry
        if optimize:
            try:
                if num_conformers > 1:
                    for conf_id in range(mol.GetNumConformers()):
                        AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
                else:
                    AllChem.MMFFOptimizeMolecule(mol)
            except Exception:
                pass  # Optimization failure is not critical
        
        # Calculate properties
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotatable = Descriptors.NumRotatableBonds(mol)
        
        # Generate output formats
        pdb_block = Chem.MolToPDBBlock(mol)
        mol_block = Chem.MolToMolBlock(mol)
        
        return {
            "smiles": smiles,
            "pdb_content": pdb_block,
            "mol_content": mol_block,
            "num_atoms": mol.GetNumAtoms(),
            "num_conformers": mol.GetNumConformers(),
            "properties": {
                "molecular_weight": round(mw, 2),
                "logP": round(logp, 2),
                "hbd": hbd,
                "hba": hba,
                "rotatable_bonds": rotatable
            }
        }
    
    @staticmethod
    def prepare_ligand_from_sdf(sdf_content: str, options: dict = None) -> Dict:
        """
        Prepare a ligand from SDF file content.
        """
        options = options or {}
        optimize = options.get('optimize', True)
        
        suppl = Chem.SDMolSupplier()
        suppl.SetData(sdf_content)
        
        mol = next(suppl, None)
        if mol is None:
            raise ValueError("Failed to parse SDF file")
        
        # Add hydrogens if not present
        mol = Chem.AddHs(mol, addCoords=True)
        
        # Optimize if requested
        if optimize:
            try:
                AllChem.MMFFOptimizeMolecule(mol)
            except Exception:
                pass
        
        smiles = Chem.MolToSmiles(mol)
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        
        return {
            "smiles": smiles,
            "pdb_content": Chem.MolToPDBBlock(mol),
            "mol_content": Chem.MolToMolBlock(mol),
            "num_atoms": mol.GetNumAtoms(),
            "properties": {
                "molecular_weight": round(mw, 2),
                "logP": round(logp, 2)
            }
        }
    
    @staticmethod
    def calculate_docking_box(pdb_content: str, ligand_residue: str = None) -> Dict:
        """
        Calculate suggested docking box from structure.
        
        If ligand_residue provided, centers on that residue.
        Otherwise, uses geometric center of protein.
        """
        from Bio.PDB import PDBParser
        import numpy as np
        
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', StringIO(pdb_content))
        
        coords = []
        ligand_coords = []
        
        for atom in structure.get_atoms():
            residue = atom.get_parent()
            res_name = residue.get_resname()
            
            if ligand_residue and res_name == ligand_residue:
                ligand_coords.append(atom.coord)
            else:
                coords.append(atom.coord)
        
        # Use ligand center if available, otherwise protein center
        if ligand_coords:
            center_coords = np.array(ligand_coords)
        else:
            center_coords = np.array(coords) if coords else np.array([[0, 0, 0]])
        
        center = center_coords.mean(axis=0)
        
        # Calculate appropriate box size
        if ligand_coords:
            # Box around ligand with padding
            ligand_coords = np.array(ligand_coords)
            extent = ligand_coords.max(axis=0) - ligand_coords.min(axis=0)
            box_size = (extent + 10).tolist()  # 5 Angstrom padding each side
        else:
            # Default box size
            box_size = [20.0, 20.0, 20.0]
        
        return {
            "center": center.tolist(),
            "size": box_size,
            "based_on": "ligand" if ligand_coords else "protein_center"
        }
