"""Builder module schemas for ligand design."""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict
import uuid


class SMILESValidateRequest(BaseModel):
    """Request to validate a SMILES string."""
    smiles: str
    
    class Config:
        json_schema_extra = {
            "example": {
                "smiles": "CC(=O)Oc1ccccc1C(=O)O"
            }
        }


class SMILESValidateResponse(BaseModel):
    """Response from SMILES validation."""
    valid: bool
    canonical_smiles: Optional[str]
    molecular_formula: Optional[str]
    molecular_weight: Optional[float]
    error: Optional[str]


class Generate3DRequest(BaseModel):
    """Request to generate 3D coordinates from SMILES."""
    smiles: str
    project_id: uuid.UUID
    method: str = Field(default="etkdg", description="Generation method")
    num_conformers: int = Field(default=10, ge=1, le=100)
    minimize: bool = Field(default=True)
    force_field: str = Field(default="mmff94", description="Force field for minimization")
    name: Optional[str] = None
    
    class Config:
        json_schema_extra = {
            "example": {
                "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "project_id": "123e4567-e89b-12d3-a456-426614174000",
                "method": "etkdg",
                "num_conformers": 10,
                "minimize": True,
                "force_field": "mmff94",
                "name": "Aspirin"
            }
        }


class ConformerInfo(BaseModel):
    """Information about a generated conformer."""
    conformer_id: int
    energy: float
    rmsd: Optional[float] = None
    file_url: str


class Generate3DResponse(BaseModel):
    """Response from 3D generation."""
    ligand_id: uuid.UUID
    smiles: str
    canonical_smiles: str
    conformers: List[ConformerInfo]
    lowest_energy_conformer_id: int
    properties: Dict


class ConformerRequest(BaseModel):
    """Request to generate multiple conformers."""
    ligand_id: uuid.UUID
    num_conformers: int = Field(default=50, ge=1, le=100)
    energy_window: float = Field(default=10.0, description="kcal/mol")
    rms_threshold: float = Field(default=0.5, description="RMSD threshold for clustering")


class ConformerResponse(BaseModel):
    """Response from conformer generation."""
    ligand_id: uuid.UUID
    conformers: List[ConformerInfo]
    num_clusters: int
