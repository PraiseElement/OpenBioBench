"""
Sequence alignment API endpoints.
"""
from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.ext.asyncio import AsyncSession
from pydantic import BaseModel
from typing import List, Optional

from app.core.database import get_db
from app.core.security import get_current_active_user
from app.models.user import User
from app.services.alignment_service import AlignmentService


router = APIRouter(prefix="/alignments", tags=["Alignments"])


class AlignmentRequest(BaseModel):
    fasta_content: str


class AlignedSequence(BaseModel):
    id: str
    sequence: str


class AlignmentResponse(BaseModel):
    aligned_sequences: List[AlignedSequence]
    num_sequences: int
    alignment_length: int
    percent_identity: float
    gap_percentage: float
    conservation: List[float]


@router.post("/align", response_model=AlignmentResponse)
async def align_sequences(
    request: AlignmentRequest,
    current_user: User = Depends(get_current_active_user),
    db: AsyncSession = Depends(get_db)
):
    """
    Perform multiple sequence alignment on FASTA sequences.
    """
    try:
        result = AlignmentService.perform_multiple_alignment(request.fasta_content)
        return result
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Alignment failed: {str(e)}"
        )
