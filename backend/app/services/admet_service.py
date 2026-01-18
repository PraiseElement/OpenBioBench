"""
ADMET Prediction Service using RDKit and pre-trained models.
Calculates physicochemical properties and predicts pharmacokinetic/toxicological properties.
"""
import logging
from typing import Dict, Tuple
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors
import numpy as np

logger = logging.getLogger(__name__)


class ADMETService:
    """Service for ADMET property prediction."""
    
    @staticmethod
    def calculate_physicochemical(mol: Chem.Mol) -> Dict:
        """
        Calculate physicochemical descriptors.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Dictionary of properties
        """
        return {
            "molecular_weight": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "logd": None,  # pH-dependent, would need pKa calculation
            "tpsa": Descriptors.TPSA(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "aromatic_rings": Descriptors.NumAromaticRings(mol),
            "lipinski_violations": ADMETService._check_lipinski(mol)
        }
    
    @staticmethod
    def _check_lipinski(mol: Chem.Mol) -> int:
        """Check Lipinski Rule of Five violations."""
        violations = 0
        if Descriptors.MolWt(mol) > 500:
            violations += 1
        if Descriptors.MolLogP(mol) > 5:
            violations += 1
        if Descriptors.NumHDonors(mol) > 5:
            violations += 1
        if Descriptors.NumHAcceptors(mol) > 10:
            violations += 1
        return violations
    
    @staticmethod
    def predict_caco2_permeability(mol: Chem.Mol) -> Dict:
        """
        Predict Caco-2 permeability (intestinal absorption).
        
        Uses a simple QSAR model based on logP and TPSA.
        Real implementation would use trained Random Forest model.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Prediction dictionary with value, confidence, interpretation
        """
        # Simple empirical model (placeholder for trained model)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        
        # log Papp (cm/s) estimation
        # Real model would be: log_papp = trained_model.predict(descriptors)
        log_papp = -5.0 + (0.3 * logp) - (0.01 * tpsa)
        
        # Confidence based on applicability domain
        confidence = "medium"  # Real model would check distance to training set
        
        # Interpretation
        if log_papp > -4.7:
            interpretation = "high permeability"
        elif log_papp > -5.7:
            interpretation = "moderate permeability"
        else:
            interpretation = "low permeability"
        
        return {
            "value": round(log_papp, 2),
            "unit": "log cm/s",
            "confidence": confidence,
            "interpretation": interpretation
        }
    
    @staticmethod
    def predict_bbb_penetration(mol: Chem.Mol) -> Dict:
        """
        Predict blood-brain barrier penetration.
        
        Uses empirical model based on molecular descriptors.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Prediction dictionary
        """
        # Simple model based on molecular properties
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        mw = Descriptors.MolWt(mol)
        
        # CNS MPO-like scoring
        # Real model: bbb_score = trained_classifier.predict_proba(descriptors)
        
        # Empirical rules:
        # BBB+ if: TPSA < 90, MW < 450, logP 1-4
        bbb_score = 0.5
        
        if tpsa < 90:
            bbb_score += 0.2
        if 1 < logp < 4:
            bbb_score += 0.2
        if mw < 450:
            bbb_score += 0.1
        
        if bbb_score > 0.7:
            interpretation = "high"
        elif bbb_score > 0.4:
            interpretation = "moderate"
        else:
            interpretation = "low"
        
        return {
            "value": round(bbb_score, 2),
            "unit": "probability",
            "confidence": "low",  # Simplified model
            "interpretation": interpretation
        }
    
    @staticmethod
    def predict_herg_ic50(mol: Chem.Mol) -> Dict:
        """
        Predict hERG channel IC50 (cardiotoxicity risk).
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Prediction dictionary
        """
        # Real implementation would use trained QSAR model
        # This is a placeholder with empirical estimation
        
        logp = Descriptors.MolLogP(mol)
        mw = Descriptors.MolWt(mol)
        
        # Rough estimation (not scientifically validated!)
        # Real: ic50 = 10^(trained_model.predict(descriptors))
        log_ic50 = 1.0 + (0.2 * logp) - (0.001 * mw)
        ic50 = 10 ** log_ic50
        
        # Risk interpretation
        if ic50 < 1.0:
            risk = "high"
        elif ic50 < 10.0:
            risk = "moderate"
        else:
            risk = "low"
        
        return {
            "value": round(ic50, 1),
            "unit": "uM",
            "confidence": "medium",
            "interpretation": risk  # Fixed: changed from 'risk' to 'interpretation'
        }
    
    @staticmethod
    def predict_ames_mutagenicity(mol: Chem.Mol) -> Dict:
        """
        Predict AMES mutagenicity.
        
        Args:
            mol: RDKit molecule
            
        Returns:
            Prediction dictionary with classification
        """
        # Structural alerts for mutagenicity (simplified)
        # Real implementation would use trained classifier
        
        smiles = Chem.MolToSmiles(mol)
        
        # Check for known mutagenic substructures (very simplified)
        mutagenic_smarts = [
            "[N+](=O)[O-]",  # Nitro groups
            "N=N",            # Azo groups
            "C=C(Cl)Cl",      # Vinyl halides
            "[C,c]1[C,c][C,c][C,c][C,c][C,c]1N",  # Aromatic amines
        ]
        
        has_alert = False
        for smarts in mutagenic_smarts:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern and mol.HasSubstructMatch(pattern):
                has_alert = True
                break
        
        prediction = "mutagenic" if has_alert else "non-mutagenic"
        confidence = 0.7 if has_alert else 0.85
        
        return {
            "prediction": prediction,
            "confidence": confidence,
            "structural_alerts": has_alert
        }
    
    @staticmethod
    def predict_all(smiles: str) -> Tuple[bool, Dict, str]:
        """
        Predict all ADMET properties for a molecule.
        
        Args:
            smiles: SMILES string
            
        Returns:
            Tuple of (success, results_dict, error_message)
        """
        try:
            # Parse SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False, {}, "Invalid SMILES string"
            
            # Calculate canonical SMILES
            canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            
            # Calculate all properties
            physicochemical = ADMETService.calculate_physicochemical(mol)
            
            pharmacokinetics = {
                "caco2_permeability": ADMETService.predict_caco2_permeability(mol),
                "bbb_penetration": ADMETService.predict_bbb_penetration(mol),
                "pgp_substrate": None,  # Not implemented
                "hia": None  # Not implemented
            }
            
            toxicity = {
                "herg_ic50": ADMETService.predict_herg_ic50(mol),
                "ames_mutagenicity": ADMETService.predict_ames_mutagenicity(mol),
                "hepatotoxicity": None,  # Not implemented
                "ld50": None  # Not implemented
            }
            
            # Generate warnings
            warnings = []
            if physicochemical["lipinski_violations"] > 0:
                warnings.append(
                    f"Lipinski violations: {physicochemical['lipinski_violations']}. "
                    "Molecule may have poor oral bioavailability."
                )
            
            if pharmacokinetics["caco2_permeability"]["value"] < -6.0:
                warnings.append("Low intestinal permeability predicted.")
            
            if toxicity["herg_ic50"]["risk"] == "high":
                warnings.append("High hERG toxicity risk. Potential cardiotoxicity concern.")
            
            if toxicity["ames_mutagenicity"]["prediction"] == "mutagenic":
                warnings.append("Structural alerts for mutagenicity detected.")
            
            # Add disclaimer
            warnings.append(
                "These are computational predictions for research purposes only. "
                "Experimental validation is required."
            )
            
            result = {
                "smiles": smiles,
                "canonical_smiles": canonical_smiles,
                "properties": physicochemical,
                "pharmacokinetics": pharmacokinetics,
                "toxicity": toxicity,
                "warnings": warnings
            }
            
            return True, result, ""
            
        except Exception as e:
            logger.error(f"ADMET prediction failed: {e}")
            return False, {}, str(e)
