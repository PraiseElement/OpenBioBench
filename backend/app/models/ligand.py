"""
Ligand model for small molecules.
Stores SMILES, structure files, and computed properties.
"""
import uuid
from datetime import datetime
from sqlalchemy import String, Text, DateTime, ForeignKey, Enum as SQLEnum
from sqlalchemy.dialects.postgresql import UUID, JSONB
from sqlalchemy.orm import Mapped, mapped_column
import enum

from app.core.database import Base


class LigandSource(str, enum.Enum):
    """Source of ligand data."""
    PUBCHEM = "pubchem"
    CHEMBL = "chembl"
    UPLOAD = "upload"
    DRAWN = "drawn"
    BUILDER = "builder"


class Ligand(Base):
    """Small molecule ligand for docking and analysis."""
    
    __tablename__ = "ligands"
    
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
    source: Mapped[LigandSource] = mapped_column(
        SQLEnum(LigandSource),
        nullable=False
    )
    external_id: Mapped[str | None] = mapped_column(String(100), nullable=True, index=True)
    
    # Chemical identifiers
    smiles: Mapped[str] = mapped_column(Text, nullable=False, index=True)
    inchi: Mapped[str | None] = mapped_column(Text, nullable=True)
    inchi_key: Mapped[str | None] = mapped_column(String(255), nullable=True, index=True)
    
    # Name and description
    name: Mapped[str | None] = mapped_column(String(255), nullable=True)
    
    # File storage (3D structure in MOL2, SDF, etc.)
    mol_file_url: Mapped[str | None] = mapped_column(String(500), nullable=True)
    
    # Computed molecular properties (MW, logP, HBD, HBA, etc.)
    properties: Mapped[dict] = mapped_column(JSONB, default=dict, nullable=False)
    
    # Additional metadata
    meta_data: Mapped[dict] = mapped_column(JSONB, default=dict, nullable=False)
    
    # Timestamps
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        default=datetime.utcnow,
        nullable=False
    )
    
    def __repr__(self) -> str:
        return f"<Ligand {self.name or self.id} ({self.smiles[:30]}...)>"
