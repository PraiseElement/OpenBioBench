"""
Ligand Builder Service using RDKit.
Handles 2D/3D molecule generation, conformer search, and energy minimization.
"""
import logging
from typing import List, Tuple, Optional, Dict
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
import numpy as np

logger = logging.getLogger(__name__)


class BuilderService:
    """Service for ligand building and 3D generation."""
    
    @staticmethod
    def validate_smiles(smiles: str) -> Tuple[bool, Optional[str], Optional[str]]:
        """
        Validate a SMILES string and return canonical form.
        
        Args:
            smiles: SMILES string to validate
            
        Returns:
            Tuple of (is_valid, canonical_smiles, error_message)
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False, None, "Invalid SMILES: could not parse molecule"
            
            # Sanitize molecule
            Chem.SanitizeMol(mol)
            
            # Generate canonical SMILES
            canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            
            return True, canonical_smiles, None
            
        except Exception as e:
            return False, None, f"SMILES validation error: {str(e)}"
    
    @staticmethod
    def calculate_properties(mol: Chem.Mol) -> Dict:
        """
        Calculate molecular properties for a molecule.
        
        Args:
            mol: RDKit molecule object
            
        Returns:
            Dictionary of calculated properties
        """
        return {
            "molecular_weight": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "tpsa": Descriptors.TPSA(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "aromatic_rings": Descriptors.NumAromaticRings(mol),
            "num_atoms": mol.GetNumAtoms(),
            "num_heavy_atoms": mol.GetNumHeavyAtoms(),
            "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
            "formal_charge": Chem.GetFormalCharge(mol),
            # Lipinski Rule of Five violations
            "lipinski_violations": sum([
                Descriptors.MolWt(mol) > 500,
                Descriptors.MolLogP(mol) > 5,
                Descriptors.NumHDonors(mol) > 5,
                Descriptors.NumHAcceptors(mol) > 10
            ])
        }
    
    @staticmethod
    def generate_3d_conformer(
        smiles: str,
        method: str = "etkdg",
        random_seed: int = 42
    ) -> Optional[Chem.Mol]:
        """
        Generate a single 3D conformer from SMILES.
        
        Args:
            smiles: SMILES string
            method: Generation method ("etkdg", "etdg", "basic")
            random_seed: Random seed for reproducibility
            
        Returns:
            RDKit molecule with 3D coordinates or None if failed
        """
        try:
            # Parse SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.error(f"Could not parse SMILES: {smiles}")
                return None
            
            # Add hydrogens
            mol = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            if method == "etkdg":
                # Experimental Torsion-angle Knowledge-based Distance Geometry
                params = AllChem.ETKDGv3()
                params.randomSeed = random_seed
                result = AllChem.EmbedMolecule(mol, params)
            elif method == "etdg":
                params = AllChem.ETDG()
                params.randomSeed = random_seed
                result = AllChem.EmbedMolecule(mol, params)
            else:
                result = AllChem.EmbedMolecule(mol, randomSeed=random_seed)
            
            if result == -1:
                logger.error(f"Could not generate 3D coordinates for: {smiles}")
                return None
            
            return mol
            
        except Exception as e:
            logger.error(f"Error generating 3D conformer: {e}")
            return None
    
    @staticmethod
    def minimize_energy(
        mol: Chem.Mol,
        force_field: str = "mmff94",
        max_iterations: int = 1000
    ) -> Tuple[Chem.Mol, float]:
        """
        Energy minimize a molecule using force field.
        
        Args:
            mol: RDKit molecule with 3D coordinates
            force_field: Force field to use ("mmff94" or "uff")
            max_iterations: Maximum optimization iterations
            
        Returns:
            Tuple of (minimized_mol, energy_kcal_mol)
        """
        try:
            if force_field.lower() == "mmff94":
                # MMFF94 force field (better for drug-like molecules)
                props = AllChem.MMFFGetMoleculeProperties(mol)
                if props is None:
                    logger.warning("MMFF94 failed, falling back to UFF")
                    force_field = "uff"
                else:
                    ff = AllChem.MMFFGetMoleculeForceField(mol, props)
                    ff.Initialize()
                    ff.Minimize(maxIts=max_iterations)
                    energy = ff.CalcEnergy()
                    return mol, energy
            
            if force_field.lower() == "uff":
                # Universal Force Field (less accurate but always works)
                ff = AllChem.UFFGetMoleculeForceField(mol)
                ff.Initialize()
                ff.Minimize(maxIts=max_iterations)
                energy = ff.CalcEnergy()
                return mol, energy
            
            raise ValueError(f"Unknown force field: {force_field}")
            
        except Exception as e:
            logger.error(f"Energy minimization failed: {e}")
            # Return molecule unchanged with NaN energy
            return mol, float('nan')
    
    @staticmethod
    def generate_conformers(
        smiles: str,
        num_conformers: int = 10,
        method: str = "etkdg",
        force_field: str = "mmff94",
        rms_threshold: float = 0.5,
        energy_window: float = 10.0,
        random_seed: int = 42
    ) -> List[Tuple[Chem.Mol, float, int]]:
        """
        Generate multiple conformers and cluster by RMSD.
        
        Args:
            smiles: SMILES string
            num_conformers: Number of conformers to generate
            method: Generation method
            force_field: Force field for minimization
            rms_threshold: RMSD threshold for clustering (Angstroms)
            energy_window: Energy window for keeping conformers (kcal/mol)
            random_seed: Random seed
            
        Returns:
            List of (mol, energy, conformer_id) tuples sorted by energy
        """
        try:
            # Parse and prepare molecule
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return []
            
            mol = Chem.AddHs(mol)
            
            # Generate conformers
            if method == "etkdg":
                params = AllChem.ETKDGv3()
                params.randomSeed = random_seed
                params.numThreads = 0  # Use all available threads
                confIds = AllChem.EmbedMultipleConfs(
                    mol,
                    numConfs=num_conformers,
                    params=params
                )
            else:
                confIds = AllChem.EmbedMultipleConfs(
                    mol,
                    numConfs=num_conformers,
                    randomSeed=random_seed
                )
            
            if len(confIds) == 0:
                logger.error(f"No conformers generated for: {smiles}")
                return []
            
            # Minimize each conformer and calculate energy
            conformers = []
            for confId in confIds:
                try:
                    if force_field.lower() == "mmff94":
                        props = AllChem.MMFFGetMoleculeProperties(mol)
                        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=confId)
                        ff.Minimize()
                        energy = ff.CalcEnergy()
                    else:
                        ff = AllChem.UFFGetMoleculeForceField(mol, confId=confId)
                        ff.Minimize()
                        energy = ff.CalcEnergy()
                    
                    conformers.append((mol, energy, confId))
                    
                except Exception as e:
                    logger.warning(f"Could not minimize conformer {confId}: {e}")
                    continue
            
            # Sort by energy
            conformers.sort(key=lambda x: x[1])
            
            # Filter by energy window
            if len(conformers) > 0:
                min_energy = conformers[0][1]
                conformers = [
                    c for c in conformers
                    if c[1] - min_energy <= energy_window
                ]
            
            # Cluster by RMSD (simple greedy clustering)
            clustered = []
            for mol_conf, energy, conf_id in conformers:
                # Check if similar to any existing cluster representative
                is_unique = True
                for existing_mol, existing_energy, existing_id in clustered:
                    # Calculate RMSD between conformers
                    rmsd = AllChem.GetConformerRMS(
                        mol,
                        conf_id,
                        existing_id,
                        prealigned=False
                    )
                    if rmsd < rms_threshold:
                        is_unique = False
                        break
                
                if is_unique:
                    clustered.append((mol_conf, energy, conf_id))
            
            logger.info(
                f"Generated {len(confIds)} conformers, "
                f"kept {len(clustered)} unique after clustering"
            )
            
            return clustered
            
        except Exception as e:
            logger.error(f"Conformer generation failed: {e}")
            return []
    
    @staticmethod
    def calculate_inchi(mol: Chem.Mol) -> Tuple[str, str]:
        """
        Calculate InChI and InChIKey for a molecule.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Tuple of (inchi, inchi_key)
        """
        try:
            inchi = Chem.MolToInchi(mol)
            inchi_key = Chem.MolToInchiKey(mol)
            return inchi, inchi_key
        except Exception as e:
            logger.error(f"InChI calculation failed: {e}")
            return "", ""
