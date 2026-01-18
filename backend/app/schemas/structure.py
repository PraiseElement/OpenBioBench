"""Structure schemas."""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict
from datetime import datetime
import uuid


class StructureCreate(BaseModel):
    """Schema for creating a structure."""
    project_id: uuid.UUID
    source: str
    pdb_id: Optional[str] = None
    file_data: Optional[str] = None  # Base64 encoded or file content
    
    class Config:
        json_schema_extra = {
            "example": {
                "project_id": "123e4567-e89b-12d3-a456-426614174000",
                "source": "pdb",
                "pdb_id": "1ABC"
            }
        }


class StructureResponse(BaseModel):
    """Schema for structure data in responses."""
    id: uuid.UUID
    project_id: uuid.UUID
    source: str
    pdb_id: Optional[str]
    file_url: str
    file_format: str
    resolution: Optional[float]
    method: str
    chains: List[str]
    num_residues: Optional[int]
    metadata: Dict
    created_at: datetime
    
    class Config:
        from_attributes = True
