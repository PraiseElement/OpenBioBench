"""ADMET prediction schemas."""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict
import uuid


class ADMETRequest(BaseModel):
    """Request for ADMET property prediction."""
    smiles: Optional[str] = None
    ligand_id: Optional[uuid.UUID] = None
    
    class Config:
        json_schema_extra = {
            "example": {
                "smiles": "CC(=O)Oc1ccccc1C(=O)O"
            }
        }


class PropertyPrediction(BaseModel):
    """Single property prediction with confidence."""
    value: float
    unit: Optional[str] = None
    confidence: str  # "low", "medium", "high"
    interpretation: Optional[str] = None


class PhysicochemicalProperties(BaseModel):
    """Physicochemical properties."""
    molecular_weight: float
    logp: float
    logd: Optional[float]
    tpsa: float
    hbd: int
    hba: int
    rotatable_bonds: int
    aromatic_rings: int
    lipinski_violations: int


class PharmacokineticProperties(BaseModel):
    """Pharmacokinetic properties."""
    caco2_permeability: PropertyPrediction
    bbb_penetration: PropertyPrediction
    pgp_substrate: Optional[PropertyPrediction]
    hia: Optional[PropertyPrediction]


class ToxicityProperties(BaseModel):
    """Toxicity predictions."""
    herg_ic50: PropertyPrediction
    ames_mutagenicity: Dict  # {"prediction": "non-mutagenic", "confidence": 0.85}
    hepatotoxicity: Optional[Dict]
    ld50: Optional[PropertyPrediction]


class ADMETResponse(BaseModel):
    """Response from ADMET prediction."""
    ligand_id: Optional[uuid.UUID]
    smiles: str
    canonical_smiles: str
    properties: PhysicochemicalProperties
    pharmacokinetics: PharmacokineticProperties
    toxicity: ToxicityProperties
    warnings: List[str] = Field(default_factory=list)
    
    class Config:
        json_schema_extra = {
            "example": {
                "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "properties": {
                    "molecular_weight": 180.16,
                    "logp": 1.19,
                    "tpsa": 63.6,
                    "hbd": 1,
                    "hba": 4,
                    "rotatable_bonds": 3,
                    "aromatic_rings": 1,
                    "lipinski_violations": 0
                },
                "pharmacokinetics": {
                    "caco2_permeability": {
                        "value": -5.2,
                        "unit": "log cm/s",
                        "confidence": "medium",
                        "interpretation": "moderate permeability"
                    },
                    "bbb_penetration": {
                        "value": 0.15,
                        "confidence": "low",
                        "interpretation": "low"
                    }
                },
                "toxicity": {
                    "herg_ic50": {
                        "value": 12.3,
                        "unit": "uM",
                        "confidence": "medium"
                    },
                    "ames_mutagenicity": {
                        "prediction": "non-mutagenic",
                        "confidence": 0.85
                    }
                },
                "warnings": []
            }
        }
