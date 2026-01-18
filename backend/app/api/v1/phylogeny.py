"""
Phylogenetic analysis API endpoints.
"""
from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.ext.asyncio import AsyncSession
from pydantic import BaseModel
from typing import List, Optional

from app.core.database import get_db
from app.core.security import get_current_active_user
from app.models.user import User
from app.services.phylogeny_service import PhylogenyService


router = APIRouter(prefix="/phylogeny", tags=["Phylogeny"])


class PhylogenyRequest(BaseModel):
    alignment_fasta: str
    method: str = "nj"  # nj (Neighbor-Joining) or upgma


class PhylogenyResponse(BaseModel):
    newick: str
    ascii_tree: Optional[str] = None
    num_taxa: int
    method: str
    total_branch_length: float


@router.post("/infer", response_model=PhylogenyResponse)
async def infer_phylogeny(
    request: PhylogenyRequest,
    current_user: User = Depends(get_current_active_user),
    db: AsyncSession = Depends(get_db)
):
    """
    Build phylogenetic tree from sequence alignment.
    Supports Neighbor-Joining and UPGMA methods.
    """
    try:
        result = PhylogenyService.build_tree(
            alignment_fasta=request.alignment_fasta,
            method=request.method
        )
        return result
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Phylogeny inference failed: {str(e)}"
        )
