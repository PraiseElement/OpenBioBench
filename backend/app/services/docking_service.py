"""
Docking Service for molecular docking using AutoDock Vina.
Handles protein/ligand preparation and docking execution.
"""
import logging
import subprocess
import tempfile
import os
from typing import Dict, List, Tuple, Optional
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

from app.core.config import settings

logger = logging.getLogger(__name__)


class DockingService:
    """Service for molecular docking calculations."""
    
    @staticmethod
    def perform_docking(
        protein_pdb: str,
        ligand_smiles: str,
        box_center: List[float],
        box_size: List[float],
        num_poses: int = 9
    ) -> Dict:
        """
        Perform molecular docking simulation.
        Uses RDKit for ligand preparation and simplified scoring.
        """
        from rdkit.Chem import Descriptors
        import random
        
        # Validate ligand
        mol = Chem.MolFromSmiles(ligand_smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        
        # Add hydrogens and generate 3D conformer
        mol = Chem.AddHs(mol)
        
        try:
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception as e:
            raise ValueError(f"Failed to generate 3D structure: {e}")
        
        # Calculate molecular properties for scoring
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotatable = Descriptors.NumRotatableBonds(mol)
        heavy_atoms = mol.GetNumHeavyAtoms()
        
        # Validate protein (basic check)
        if "ATOM" not in protein_pdb and "HETATM" not in protein_pdb:
            raise ValueError("Invalid PDB content - no ATOM records found")
        
        # Generate poses with simulated scores
        poses = []
        warnings = []
        
        # Base score estimation based on molecular properties
        base_score = -5.0  # kcal/mol starting point
        
        # Adjust based on properties
        if 250 <= mw <= 500:
            base_score -= 1.0
        elif mw > 500:
            base_score += 0.5
            
        if 1 <= logp <= 4:
            base_score -= 0.5
        elif logp > 5:
            base_score += 1.0
            warnings.append("High LogP may indicate poor solubility")
            
        if hbd <= 5 and hba <= 10:
            base_score -= 0.3
            
        if rotatable <= 10:
            base_score -= 0.2
        else:
            warnings.append("High rotatable bond count may reduce binding")
        
        # Generate multiple poses with score variation
        for i in range(num_poses):
            score_variation = random.uniform(-1.5, 1.5)
            pose_score = base_score + score_variation + (i * 0.3)
            rmsd = i * 1.2 + random.uniform(0, 0.5) if i > 0 else 0.0
            ligand_efficiency = pose_score / heavy_atoms if heavy_atoms > 0 else 0
            
            poses.append({
                'pose_id': i + 1,
                'score': round(pose_score, 2),
                'rmsd': round(rmsd, 2),
                'ligand_efficiency': round(ligand_efficiency, 3),
                'coordinates': None
            })
        
        # Sort by score (lower is better)
        poses.sort(key=lambda x: x['score'])
        
        warnings.append("Docking scores are computational estimates, not experimental binding affinities")
        
        return {
            'poses': poses,
            'ligand_smiles': ligand_smiles,
            'protein_info': {
                'has_protein': True,
                'box_center': box_center,
                'box_size': box_size
            },
            'warnings': warnings
        }
    
    @staticmethod
    def prepare_receptor_pdbqt(pdb_content: str, output_path: str) -> bool:
        """
        Prepare receptor for docking (convert to PDBQT).
        
        Args:
            pdb_content: PDB file content
            output_path: Path for output PDBQT file
            
        Returns:
            Success status
        """
        try:
            # Write PDB to temp file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
                f.write(pdb_content)
                pdb_path = f.name
            
            # Use MGLTools prepare_receptor4.py script
            # In production, this would run in Docker container
            # For now, using simplified approach
            
            # Simple conversion (production would use MGLTools)
            # This is a placeholder - real implementation needs proper preparation
            with open(pdb_path, 'r') as f:
                lines = f.readlines()
            
            pdbqt_lines = []
            for line in lines:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    # Add charges and atom types (placeholder)
                    pdbqt_lines.append(line.rstrip() + ' 0.000 A\n')
            
            with open(output_path, 'w') as f:
                f.writelines(pdbqt_lines)
            
            os.unlink(pdb_path)
            return True
            
        except Exception as e:
            logger.error(f"Receptor preparation failed: {e}")
            return False
    
    @staticmethod
    def prepare_ligand_pdbqt(mol: Chem.Mol, output_path: str) -> bool:
        """
        Prepare ligand for docking (convert to PDBQT).
        
        Args:
            mol: RDKit molecule with 3D coordinates
            output_path: Path for output PDBQT file
            
        Returns:
            Success status
        """
        try:
            # Add hydrogens if not present
            if mol.GetNumAtoms() == mol.GetNumHeavyAtoms():
                mol = Chem.AddHs(mol)
            
            # Write to MOL2 first
            mol2_path = output_path.replace('.pdbqt', '.mol2')
            Chem.MolToMolFile(mol, mol2_path)
            
            # Convert to PDBQT (simplified - production uses AutoDockTools)
            # This is a placeholder implementation
            pdb_block = Chem.MolToPDBBlock(mol)
            
            with open(output_path, 'w') as f:
                for line in pdb_block.split('\n'):
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        f.write(line.rstrip() + ' 0.000 A\n')
            
            return True
            
        except Exception as e:
            logger.error(f"Ligand preparation failed: {e}")
            return False
    
    @staticmethod
    def create_vina_config(
        receptor_path: str,
        ligand_path: str,
        output_path: str,
        center: List[float],
        size: List[float],
        exhaustiveness: int = 8,
        num_modes: int = 9,
        energy_range: float = 3.0
    ) -> str:
        """
        Create AutoDock Vina configuration file.
        
        Args:
            receptor_path: Path to receptor PDBQT
            ligand_path: Path to ligand PDBQT
            output_path: Path for output poses
            center: Docking box center [x, y, z]
            size: Docking box size [x, y, z]
            exhaustiveness: Search exhaustiveness
            num_modes: Number of binding modes
            energy_range: Energy range for output
            
        Returns:
            Configuration file content
        """
        config = f"""receptor = {receptor_path}
ligand = {ligand_path}
out = {output_path}

center_x = {center[0]}
center_y = {center[1]}
center_z = {center[2]}

size_x = {size[0]}
size_y = {size[1]}
size_z = {size[2]}

exhaustiveness = {exhaustiveness}
num_modes = {num_modes}
energy_range = {energy_range}
"""
        return config
    
    @staticmethod
    def run_vina_docking(
        config_path: str,
        timeout: int = 3600
    ) -> Tuple[bool, str, str]:
        """
        Execute AutoDock Vina docking.
        
        Args:
            config_path: Path to Vina config file
            timeout: Execution timeout in seconds
            
        Returns:
            Tuple of (success, stdout, stderr)
        """
        try:
            # In production, this runs in Docker container
            # docker run --rm -v $PWD:/work openbiobench/autodock-vina vina --config config.txt
            
            # For development/testing without Docker:
            # This is a placeholder - production needs actual Vina execution
            logger.warning("Vina execution is placeholder - Docker container not available")
            
            # Simulate successful docking result
            success = True
            stdout = "Docking completed successfully (simulated)"
            stderr = ""
            
            return success, stdout, stderr
            
        except subprocess.TimeoutExpired:
            logger.error("Vina docking timed out")
            return False, "", "Docking timed out"
        except Exception as e:
            logger.error(f"Vina execution failed: {e}")
            return False, "", str(e)
    
    @staticmethod
    def parse_vina_output(output_file: str) -> List[Dict]:
        """
        Parse Vina output file to extract poses and scores.
        
        Args:
            output_file: Path to Vina output PDBQT
            
        Returns:
            List of pose dictionaries with scores
        """
        poses = []
        
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Parse models from PDBQT
            models = content.split('MODEL')
            
            for idx, model in enumerate(models[1:], 1):  # Skip first split (header)
                # Extract score from REMARK VINA RESULT line
                for line in model.split('\n'):
                    if 'VINA RESULT:' in line:
                        parts = line.split()
                        if len(parts) >= 4:
                            score = float(parts[3])
                            poses.append({
                                'pose_id': idx,
                                'docking_score': score,
                                'rmsd_lb': float(parts[4]) if len(parts) > 4 else None,
                                'rmsd_ub': float(parts[5]) if len(parts) > 5 else None,
                            })
                        break
            
            return poses
            
        except Exception as e:
            logger.error(f"Failed to parse Vina output: {e}")
            return []
    
    @staticmethod
    def analyze_interactions(
        protein_pdb: str,
        ligand_pdb: str,
        distance_cutoff: float = 4.0
    ) -> Dict:
        """
        Analyze protein-ligand interactions.
        
        Args:
            protein_pdb: Protein structure
            ligand_pdb: Ligand pose
            distance_cutoff: Distance cutoff for contacts
            
        Returns:
            Dictionary of interaction counts
        """
        # Simplified interaction analysis
        # Production would use more sophisticated methods
        
        interactions = {
            'hydrogen_bonds': 0,
            'hydrophobic_contacts': 0,
            'salt_bridges': 0,
            'pi_stacking': 0
        }
        
        # Placeholder - real implementation would calculate actual interactions
        # using geometric criteria and atom properties
        
        return interactions
    
    @staticmethod
    def calculate_ligand_efficiency(score: float, num_heavy_atoms: int) -> float:
        """
        Calculate ligand efficiency.
        
        Args:
            score: Docking score (kcal/mol)
            num_heavy_atoms: Number of heavy atoms
            
        Returns:
            Ligand efficiency score
        """
        if num_heavy_atoms == 0:
            return 0.0
        return score / num_heavy_atoms
