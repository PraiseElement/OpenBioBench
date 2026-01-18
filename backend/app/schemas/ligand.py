"""Ligand schemas."""
from pydantic import BaseModel, Field
from typing import Optional, Dict
from datetime import datetime
import uuid


class LigandCreate(BaseModel):
    """Schema for creating a ligand."""
    project_id: uuid.UUID
    smiles: str = Field(..., min_length=1)
    name: Optional[str] = None
    source: str = "drawn"
    
    class Config:
        json_schema_extra = {
            "example": {
                "project_id": "123e4567-e89b-12d3-a456-426614174000",
                "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "name": "Aspirin",
                "source": "drawn"
            }
        }


class LigandResponse(BaseModel):
    """Schema for ligand data in responses."""
    id: uuid.UUID
    project_id: uuid.UUID
    source: str
    smiles: str
    inchi: Optional[str]
    inchi_key: Optional[str]
    name: Optional[str]
    mol_file_url: Optional[str]
    properties: Dict
    metadata: Dict
    created_at: datetime
    
    class Config:
        from_attributes = True
