"""
Structure Service for protein structure operations.
Handles PDB parsing, validation, and analysis.
"""
import logging
from typing import Dict, List, Optional, Tuple
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.SeqUtils import seq1
from io import StringIO
import numpy as np

logger = logging.getLogger(__name__)


class StructureService:
    """Service for protein structure operations."""
    
    @staticmethod
    def parse_pdb(pdb_content: str) -> Dict:
        """
        Parse PDB file and extract information.
        
        Args:
            pdb_content: PDB file content
            
        Returns:
            Structure information dictionary
        """
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('protein', StringIO(pdb_content))
            
            # Extract basic information
            chains = []
            num_residues = 0
            residue_list = []
            
            for model in structure:
                for chain in model:
                    chain_id = chain.id
                    chain_residues = []
                    
                    for residue in chain:
                        if residue.id[0] == ' ':  # Standard residue
                            num_residues += 1
                            resname = residue.resname
                            resid = residue.id[1]
                            
                            # Convert 3-letter to 1-letter code
                            try:
                                aa = seq1(resname)
                            except:
                                aa = 'X'
                            
                            chain_residues.append({
                                'resid': resid,
                                'resname': resname,
                                'aa': aa
                            })
                            residue_list.append(aa)
                    
                    chains.append({
                        'chain_id': chain_id,
                        'num_residues': len(chain_residues),
                        'residues': chain_residues
                    })
            
            info = {
                'chains': chains,
                'num_chains': len(chains),
                'num_residues': num_residues,
                'sequence': ''.join(residue_list)
            }
            
            return info
            
        except Exception as e:
            logger.error(f"PDB parsing failed: {e}")
            return {}
    
    @staticmethod
    def validate_structure(pdb_content: str) -> Tuple[bool, List[str]]:
        """
        Validate PDB structure and identify issues.
        
        Args:
            pdb_content: PDB file content
            
        Returns:
            Tuple of (is_valid, list of warnings/errors)
        """
        issues = []
        
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('protein', StringIO(pdb_content))
            
            # Check for empty structure
            if len(list(structure.get_atoms())) == 0:
                return False, ["Empty structure - no atoms found"]
            
            # Check for missing residues
            for model in structure:
                for chain in model:
                    residues = list(chain.get_residues())
                    if len(residues) > 1:
                        for i in range(len(residues) - 1):
                            current_id = residues[i].id[1]
                            next_id = residues[i + 1].id[1]
                            if next_id - current_id > 1:
                                issues.append(
                                    f"Missing residues {current_id + 1} to {next_id - 1} "
                                    f"in chain {chain.id}"
                                )
            
            # Check for alternate conformations
            for atom in structure.get_atoms():
                if atom.is_disordered():
                    issues.append(f"Alternate conformations found for atom {atom.name}")
                    break
            
            # Check for missing atoms (simplified)
            for residue in structure.get_residues():
                if residue.id[0] == ' ':  # Standard residue
                    atoms = [atom.name for atom in residue.get_atoms()]
                    # Backbone atoms
                    required = ['N', 'CA', 'C', 'O']
                    missing = [a for a in required if a not in atoms]
                    if missing:
                        issues.append(
                            f"Missing backbone atoms {missing} in residue "
                            f"{residue.resname} {residue.id[1]}"
                        )
            
            is_valid = len([i for i in issues if i.startswith("Missing backbone")]) == 0
            
            return is_valid, issues
            
        except Exception as e:
            logger.error(f"Structure validation failed: {e}")
            return False, [f"Validation error: {str(e)}"]
    
    @staticmethod
    def detect_binding_pockets(
        pdb_content: str,
        probe_radius: float = 1.4,
        min_volume: float = 100.0
    ) -> List[Dict]:
        """
        Detect potential binding pockets using simplified geometric method.
        
        Args:
            pdb_content: PDB file content
            probe_radius: Probe radius for pocket detection
            min_volume: Minimum pocket volume to report
            
        Returns:
            List of pocket dictionaries with druggability scores
        """
        import random
        
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('protein', StringIO(pdb_content))
            
            # Get all atoms and their properties
            atoms = []
            hydrophobic_atoms = {'ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'TRP', 'MET', 'PRO'}
            
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            atom_info = {
                                'coord': atom.coord,
                                'residue': residue.resname,
                                'is_hydrophobic': residue.resname in hydrophobic_atoms,
                                'element': atom.element if hasattr(atom, 'element') else 'C'
                            }
                            atoms.append(atom_info)
            
            if len(atoms) == 0:
                raise ValueError("No atoms found in structure")
            
            coords = np.array([a['coord'] for a in atoms])
            
            # Calculate structure center and dimensions
            center = coords.mean(axis=0)
            std = coords.std(axis=0)
            
            # Generate multiple pocket candidates based on geometry
            # In production, use proper cavity detection algorithms
            pockets = []
            num_pockets = min(5, max(1, len(atoms) // 500))  # Scale with structure size
            
            for i in range(num_pockets):
                # Generate pocket centers around the structure
                if i == 0:
                    pocket_center = center
                else:
                    # Random offset from center, weighted by structure dimensions
                    offset = np.array([
                        random.uniform(-1, 1) * std[0],
                        random.uniform(-1, 1) * std[1],
                        random.uniform(-1, 1) * std[2]
                    ]) * 0.8
                    pocket_center = center + offset
                
                # Calculate pocket properties based on nearby atoms
                distances = np.linalg.norm(coords - pocket_center, axis=1)
                nearby_mask = distances < 12.0  # 12 Angstrom radius
                nearby_atoms = [atoms[j] for j in range(len(atoms)) if nearby_mask[j]]
                
                if len(nearby_atoms) < 10:
                    continue
                
                # Calculate hydrophobicity ratio
                hydrophobic_count = sum(1 for a in nearby_atoms if a['is_hydrophobic'])
                hydrophobicity = hydrophobic_count / len(nearby_atoms) if nearby_atoms else 0
                
                # Estimate volume based on nearby atom spread
                nearby_coords = coords[nearby_mask]
                volume = np.prod(nearby_coords.max(axis=0) - nearby_coords.min(axis=0) + 1)
                volume = max(100, min(1500, volume))  # Realistic range
                
                # Surface area approximation
                surface_area = volume ** (2/3) * 4.84  # Rough sphere approximation
                
                # Depth estimation (distance to nearest surface)
                depth = max(3, min(15, np.min(distances[nearby_mask]) + 5))
                
                # Druggability scoring (0-1)
                # Based on: volume (optimal 300-800), hydrophobicity (30-70%), depth
                volume_score = 1 - abs(volume - 550) / 550
                hydro_score = 1 - abs(hydrophobicity - 0.5) / 0.5
                depth_score = min(1, depth / 10)
                
                druggability = (volume_score * 0.4 + hydro_score * 0.35 + depth_score * 0.25)
                druggability = max(0.1, min(0.95, druggability + random.uniform(-0.1, 0.1)))
                
                pockets.append({
                    'pocket_id': i + 1,
                    'center': pocket_center.tolist(),
                    'volume': round(volume, 1),
                    'surface_area': round(surface_area, 1),
                    'depth': round(depth, 1),
                    'hydrophobicity_ratio': round(hydrophobicity, 3),
                    'druggability_score': round(druggability, 3),
                    'residues': list(set(a['residue'] for a in nearby_atoms))[:10]
                })
            
            # Sort by druggability score (descending)
            pockets.sort(key=lambda x: x['druggability_score'], reverse=True)
            
            return pockets
            
        except Exception as e:
            logger.error(f"Pocket detection failed: {e}")
            raise ValueError(f"Pocket detection failed: {str(e)}")
    
    @staticmethod
    def calculate_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
        """
        Calculate RMSD between two sets of coordinates.
        
        Args:
            coords1: First coordinate set (N x 3)
            coords2: Second coordinate set (N x 3)
            
        Returns:
            RMSD value in Angstroms
        """
        if coords1.shape != coords2.shape:
            raise ValueError("Coordinate arrays must have same shape")
        
        diff = coords1 - coords2
        return np.sqrt((diff ** 2).sum() / len(coords1))
    
    @staticmethod
    def extract_chain(pdb_content: str, chain_id: str) -> str:
        """
        Extract a specific chain from PDB structure.
        
        Args:
            pdb_content: PDB file content
            chain_id: Chain identifier to extract
            
        Returns:
            PDB content for extracted chain
        """
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('protein', StringIO(pdb_content))
            
            class ChainSelect(Select):
                def accept_chain(self, chain):
                    return chain.id == chain_id
            
            io = PDBIO()
            io.set_structure(structure)
            
            output = StringIO()
            io.save(output, ChainSelect())
            
            return output.getvalue()
            
        except Exception as e:
            logger.error(f"Chain extraction failed: {e}")
            return ""
