"""Docking module schemas."""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict
import uuid


class DockingBox(BaseModel):
    """Docking box configuration."""
    center: List[float] = Field(..., min_length=3, max_length=3, description="X, Y, Z coordinates")
    size: List[float] = Field(..., min_length=3, max_length=3, description="Box size in Angstroms")
    
    class Config:
        json_schema_extra = {
            "example": {
                "center": [12.5, -3.2, 45.8],
                "size": [20.0, 20.0, 20.0]
            }
        }


class DockingRequest(BaseModel):
    """Request to run molecular docking."""
    project_id: uuid.UUID
    protein_structure_id: uuid.UUID
    ligand_ids: List[uuid.UUID] = Field(..., min_length=1)
    docking_box: DockingBox
    exhaustiveness: int = Field(default=8, ge=1, le=32)
    num_poses: int = Field(default=9, ge=1, le=20)
    energy_range: float = Field(default=3.0, ge=1.0, le=10.0)
    
    class Config:
        json_schema_extra = {
            "example": {
                "project_id": "123e4567-e89b-12d3-a456-426614174000",
                "protein_structure_id": "789e4567-e89b-12d3-a456-426614174111",
                "ligand_ids": ["456e4567-e89b-12d3-a456-426614174222"],
                "docking_box": {
                    "center": [12.5, -3.2, 45.8],
                    "size": [20.0, 20.0, 20.0]
                },
                "exhaustiveness": 8,
                "num_poses": 9,
                "energy_range": 3.0
            }
        }


class InteractionInfo(BaseModel):
    """Protein-ligand interaction information."""
    hydrogen_bonds: int
    hydrophobic_contacts: int
    salt_bridges: int
    pi_stacking: int


class DockingPose(BaseModel):
    """Single docking pose result."""
    pose_id: int
    docking_score: float
    ligand_efficiency: float
    file_url: str
    interactions: InteractionInfo


class DockingResult(BaseModel):
    """Docking result for one ligand."""
    ligand_id: uuid.UUID
    ligand_name: Optional[str]
    poses: List[DockingPose]
    best_score: float


class DockingResponse(BaseModel):
    """Response from docking job."""
    job_id: uuid.UUID
    status: str
    completed_at: Optional[str]
    results: List[DockingResult]
    warnings: List[str] = Field(default_factory=list)
    provenance: Dict
