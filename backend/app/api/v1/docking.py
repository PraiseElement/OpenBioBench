"""
Molecular docking API endpoints.
"""
from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.ext.asyncio import AsyncSession
from pydantic import BaseModel
from typing import List, Optional

from app.core.database import get_db
from app.core.security import get_current_active_user
from app.models.user import User
from app.services.docking_service import DockingService


router = APIRouter(prefix="/docking", tags=["Docking"])


class DockingRequest(BaseModel):
    protein_pdb: str
    ligand_smiles: str
    box_center: List[float] = [0, 0, 0]
    box_size: List[float] = [20, 20, 20]
    exhaustiveness: int = 8
    num_poses: int = 9


class DockingPose(BaseModel):
    pose_id: int
    score: float
    rmsd: float
    ligand_efficiency: float
    coordinates: Optional[List[List[float]]] = None


class DockingResponse(BaseModel):
    poses: List[DockingPose]
    ligand_smiles: str
    protein_info: Optional[dict] = None
    warnings: List[str] = []


@router.post("/run", response_model=DockingResponse)
async def run_docking(
    request: DockingRequest,
    current_user: User = Depends(get_current_active_user),
    db: AsyncSession = Depends(get_db)
):
    """
    Run molecular docking simulation.
    Returns ranked binding poses with scores.
    """
    try:
        result = DockingService.perform_docking(
            protein_pdb=request.protein_pdb,
            ligand_smiles=request.ligand_smiles,
            box_center=request.box_center,
            box_size=request.box_size,
            num_poses=request.num_poses
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
            detail=f"Docking failed: {str(e)}"
        )
