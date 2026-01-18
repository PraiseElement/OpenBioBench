"""
Binding pocket detection API endpoints.
"""
from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.ext.asyncio import AsyncSession
from pydantic import BaseModel
from typing import List, Optional

from app.core.database import get_db
from app.core.security import get_current_active_user
from app.models.user import User
from app.services.structure_service import StructureService


router = APIRouter(prefix="/pockets", tags=["Pockets"])


class PocketDetectionRequest(BaseModel):
    pdb_content: str


class Pocket(BaseModel):
    pocket_id: int
    center: List[float]
    volume: float
    surface_area: float
    depth: float
    hydrophobicity_ratio: float
    druggability_score: float
    residues: Optional[List[str]] = None


class PocketDetectionResponse(BaseModel):
    pockets: List[Pocket]
    structure_info: Optional[dict] = None


@router.post("/detect", response_model=PocketDetectionResponse)
async def detect_pockets(
    request: PocketDetectionRequest,
    current_user: User = Depends(get_current_active_user),
    db: AsyncSession = Depends(get_db)
):
    """
    Detect binding pockets in a protein structure.
    Returns ranked pockets with druggability scores.
    """
    try:
        pockets = StructureService.detect_binding_pockets(request.pdb_content)
        return {"pockets": pockets}
    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Pocket detection failed: {str(e)}"
        )
