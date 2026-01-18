"""
Result model for storing computational outputs.
Links to jobs and stores file references with summaries.
"""
import uuid
from datetime import datetime
from sqlalchemy import String, DateTime, ForeignKey, Text
from sqlalchemy.dialects.postgresql import UUID, JSONB
from sqlalchemy.orm import Mapped, mapped_column

from app.core.database import Base


class Result(Base):
    """Computational result from a job."""
    
    __tablename__ = "results"
    
    id: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True),
        primary_key=True,
        default=uuid.uuid4
    )
    job_id: Mapped[uuid.UUID] = mapped_column(
        UUID(as_uuid=True),
        ForeignKey("jobs.id", ondelete="CASCADE"),
        nullable=False,
        index=True
    )
    
    # Result type (e.g., "alignment", "docking_pose", "phylo_tree", "admet_report")
    result_type: Mapped[str] = mapped_column(String(100), nullable=False)
    
    # File storage (S3/MinIO path)
    file_url: Mapped[str | None] = mapped_column(String(500), nullable=True)
    
    # Quick summary for display (key metrics, scores, etc.)
    summary: Mapped[dict] = mapped_column(JSONB, default=dict, nullable=False)
    
    # Full data (for complex results that don't need separate files)
    data: Mapped[dict | None] = mapped_column(JSONB, nullable=True)
    
    # Timestamps
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        default=datetime.utcnow,
        nullable=False
    )
    
    def __repr__(self) -> str:
        return f"<Result {self.result_type} for Job {self.job_id}>"
