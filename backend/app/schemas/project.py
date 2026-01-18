"""Project schemas."""
from pydantic import BaseModel, Field
from typing import Optional, List
from datetime import datetime
import uuid


class ProjectCreate(BaseModel):
    """Schema for creating a new project."""
    name: str = Field(..., min_length=1, max_length=255)
    description: Optional[str] = None
    tags: List[str] = Field(default_factory=list)
    
    class Config:
        json_schema_extra = {
            "example": {
                "name": "COVID-19 Drug Discovery",
                "description": "Virtual screening against SARS-CoV-2 main protease",
                "tags": ["covid", "antiviral", "docking"]
            }
        }


class ProjectUpdate(BaseModel):
    """Schema for updating a project."""
    name: Optional[str] = Field(None, min_length=1, max_length=255)
    description: Optional[str] = None
    tags: Optional[List[str]] = None


class ProjectResponse(BaseModel):
    """Schema for project data in responses."""
    id: uuid.UUID
    owner_id: uuid.UUID
    name: str
    description: Optional[str]
    tags: List[str]
    created_at: datetime
    updated_at: datetime
    
    class Config:
        from_attributes = True
