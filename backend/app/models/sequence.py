"""
Sequence model for biological sequences (protein, DNA, RNA).
Supports import from UniProt, RefSeq, and manual upload.
"""
import uuid
from datetime import datetime
from sqlalchemy import String, Text, DateTime, ForeignKey, Enum as SQLEnum
from sqlalchemy.dialects.postgresql import UUID, JSONB
from sqlalchemy.orm import Mapped, mapped_column
import enum

from app.core.database import Base


class SequenceSource(str, enum.Enum):
    """Source of sequence data."""
    UPLOAD = "upload"
    UNIPROT = "uniprot"
    REFSEQ = "refseq"
    SWISS_PROT = "swiss_prot"
    EMBL = "embl"


class Sequence(Base):
    """Biological sequence (protein or nucleotide)."""
    
    __tablename__ = "sequences"
    
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
    source: Mapped[SequenceSource] = mapped_column(
        SQLEnum(SequenceSource),
        nullable=False
    )
    accession: Mapped[str | None] = mapped_column(String(100), nullable=True, index=True)
    
    # Sequence data
    sequence: Mapped[str] = mapped_column(Text, nullable=False)
    sequence_type: Mapped[str] = mapped_column(String(20), default="protein", nullable=False)
    
    # Metadata (organism, gene name, annotations, etc.)
    meta_data: Mapped[dict] = mapped_column(JSONB, default=dict, nullable=False)
    
    # Timestamps
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        default=datetime.utcnow,
        nullable=False
    )
    
    def __repr__(self) -> str:
        return f"<Sequence {self.accession or self.id} ({len(self.sequence)} residues)>"
