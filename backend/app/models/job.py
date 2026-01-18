"""
Job model for tracking computational tasks.
Implements comprehensive provenance tracking as per SDD requirements.
"""
import uuid
from datetime import datetime
from sqlalchemy import String, DateTime, ForeignKey, Integer, Enum as SQLEnum
from sqlalchemy.dialects.postgresql import UUID, JSONB, ARRAY
from sqlalchemy.orm import Mapped, mapped_column
import enum
from typing import List

from app.core.database import Base


class JobType(str, enum.Enum):
    """Type of computational job."""
    # Sequence analysis
    ALIGNMENT = "alignment"
    PHYLOGENY = "phylogeny"
    CONSERVATION = "conservation"
    
    # Structure analysis
    POCKET_DETECTION = "pocket_detection"
    STRUCTURE_PREPARATION = "structure_preparation"
    
    # Molecular operations
    LIGAND_BUILDER_2D = "ligand_builder_2d"
    LIGAND_BUILDER_3D = "ligand_builder_3d"
    LIGAND_PREPARATION = "ligand_preparation"
    
    # Docking and simulation
    DOCKING = "docking"
    MD_SIMULATION = "md_simulation"
    MD_ANALYSIS = "md_analysis"
    
    # Prediction
    ADMET_PREDICTION = "admet_prediction"
    DECISION_SUPPORT = "decision_support"
    
    # Workflow
    WORKFLOW = "workflow"


class JobStatus(str, enum.Enum):
    """Status of computational job."""
    PENDING = "pending"
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class Job(Base):
    """Computational job with full provenance tracking."""
    
    __tablename__ = "jobs"
    
    id: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True),
        primary_key=True,
        default=uuid.uuid4
    )
    user_id: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True),
        ForeignKey("users.id", ondelete="CASCADE"),
        nullable=False,
        index=True
    )
    project_id: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True),
        ForeignKey("projects.id", ondelete="CASCADE"),
        nullable=False,
        index=True
    )
    
    # Job identification
    job_type: Mapped[JobType] = mapped_column(
        SQLEnum(JobType),
        nullable=False,
        index=True
    )
    status: Mapped[JobStatus] = mapped_column(
        SQLEnum(JobStatus),
        nullable=False,
        default=JobStatus.PENDING,
        index=True
    )
    
    # Job configuration
    parameters: Mapped[dict] = mapped_column(JSONB, default=dict, nullable=False)
    
    # Input/Output tracking (UUIDs of sequences, structures, ligands)
    inputs: Mapped[List[str]] = mapped_column(ARRAY(String), default=list, nullable=False)
    outputs: Mapped[List[str]] = mapped_column(ARRAY(String), default=list, nullable=False)
    
    # Execution tracking
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        default=datetime.utcnow,
        nullable=False,
        index=True
    )
    started_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    completed_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    
    # Resource usage
    compute_time_seconds: Mapped[int | None] = mapped_column(Integer, nullable=True)
    
    # Error tracking
    error_message: Mapped[str | None] = mapped_column(String(1000), nullable=True)
    
    # Provenance (tool version, container image, command, etc.)
    provenance: Mapped[dict] = mapped_column(JSONB, default=dict, nullable=False)
    
    # Progress tracking (0-100)
    progress: Mapped[int] = mapped_column(Integer, default=0, nullable=False)
    
    def __repr__(self) -> str:
        return f"<Job {self.job_type} {self.status}>"
