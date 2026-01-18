"""
Structure model for 3D protein structures.
Supports PDB, AlphaFold, and user uploads.
"""
import uuid
from datetime import datetime
from sqlalchemy import String, Text, DateTime, ForeignKey, Float, Enum as SQLEnum
from sqlalchemy.dialects.postgresql import UUID, JSONB, ARRAY
from sqlalchemy.orm import Mapped, mapped_column
import enum
from typing import List

from app.core.database import Base


class StructureSource(str, enum.Enum):
    """Source of structure data."""
    PDB = "pdb"
    ALPHAFOLD = "alphafold"
    UPLOAD = "upload"
    HOMOLOGY_MODEL = "homology_model"


class StructureMethod(str, enum.Enum):
    """Experimental method used to determine structure."""
    XRAY = "xray"
    NMR = "nmr"
    EM = "em"
    PREDICTION = "prediction"


class Structure(Base):
    """3D protein/molecular structure."""
    
    __tablename__ = "structures"
    
    id: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True),
        primary_key=True,
        default=uuid.uuid4
    )
    project_id: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True),
        ForeignKey("projects.id", ondelete="CASCADE"),
        nullable=False,
        index=True
    )
    
    # Source information
    source: Mapped[StructureSource] = mapped_column(
        SQLEnum(StructureSource),
        nullable=False
    )
    pdb_id: Mapped[str | None] = mapped_column(String(10), nullable=True, index=True)
    
    # File storage (S3/MinIO path)
    file_url: Mapped[str] = mapped_column(String(500), nullable=False)
    file_format: Mapped[str] = mapped_column(String(20), default="pdb", nullable=False)
    
    # Quality metrics
    resolution: Mapped[float | None] = mapped_column(Float, nullable=True)
    method: Mapped[StructureMethod] = mapped_column(
        SQLEnum(StructureMethod),
        nullable=False
    )
    
    # Structural information
    chains: Mapped[List[str]] = mapped_column(ARRAY(String), default=list, nullable=False)
    num_residues: Mapped[int | None] = mapped_column(nullable=True)
    
    # Additional metadata (organism, ligands present, etc.)
    meta_data: Mapped[dict] = mapped_column(JSONB, default=dict, nullable=False)
    
    # Timestamps
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        default=datetime.utcnow,
        nullable=False
    )
    
    def __repr__(self) -> str:
        return f"<Structure {self.pdb_id or self.id} ({self.method})>"
