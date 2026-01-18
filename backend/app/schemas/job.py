"""Job schemas."""
from pydantic import BaseModel, Field
from typing import Optional, List, Dict
from datetime import datetime
import uuid


class JobCreate(BaseModel):
    """Schema for creating a job."""
    project_id: uuid.UUID
    job_type: str
    parameters: Dict = Field(default_factory=dict)
    inputs: List[str] = Field(default_factory=list)
    
    class Config:
        json_schema_extra = {
            "example": {
                "project_id": "123e4567-e89b-12d3-a456-426614174000",
                "job_type": "docking",
                "parameters": {
                    "exhaustiveness": 8,
                    "num_poses": 9
                },
                "inputs": ["protein_id", "ligand_id"]
            }
        }


class JobUpdate(BaseModel):
    """Schema for updating job status."""
    status: Optional[str] = None
    progress: Optional[int] = Field(None, ge=0, le=100)
    error_message: Optional[str] = None


class JobResponse(BaseModel):
    """Schema for job data in responses."""
    id: uuid.UUID
    user_id: uuid.UUID
    project_id: uuid.UUID
    job_type: str
    status: str
    parameters: Dict
    inputs: List[str]
    outputs: List[str]
    created_at: datetime
    started_at: Optional[datetime]
    completed_at: Optional[datetime]
    compute_time_seconds: Optional[int]
    error_message: Optional[str]
    provenance: Dict
    progress: int
    
    class Config:
        from_attributes = True
