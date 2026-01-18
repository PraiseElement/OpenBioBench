"""
ADMET prediction API endpoints.
"""
from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy import select

from app.core.database import get_db
from app.core.security import get_current_active_user
from app.models.user import User
from app.models.ligand import Ligand
from app.services.admet_service import ADMETService
from app.schemas.admet import ADMETRequest, ADMETResponse

router = APIRouter(prefix="/admet", tags=["ADMET Prediction"])


@router.post("/predict", response_model=ADMETResponse)
async def predict_admet(
    request: ADMETRequest,
    current_user: User = Depends(get_current_active_user),
    db: AsyncSession = Depends(get_db)
):
    """
    Predict ADMET properties for a molecule.
    
    Provide either a SMILES string or a ligand ID from the database.
    
    Returns:
    - Physicochemical properties (MW, logP, TPSA, etc.)
    - Pharmacokinetic predictions (Caco-2, BBB, etc.)
    - Toxicity predictions (hERG, AMES, etc.)
    - Warnings and disclaimers
    """
    smiles = request.smiles
    
    # If ligand_id provided, fetch from database
    if request.ligand_id and not smiles:
        result = await db.execute(
            select(Ligand).where(Ligand.id == request.ligand_id)
        )
        ligand = result.scalar_one_or_none()
        
        if not ligand:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Ligand not found"
            )
        
        smiles = ligand.smiles
    
    if not smiles:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Either SMILES or ligand_id must be provided"
        )
    
    # Predict ADMET properties
    success, result, error = ADMETService.predict_all(smiles)
    
    if not success:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=error
        )
    
    # Convert to response schema
    from app.schemas.admet import (
        PhysicochemicalProperties,
        PharmacokineticProperties,
        ToxicityProperties
    )
    
    admet_response = ADMETResponse(
        ligand_id=request.ligand_id,
        smiles=result["smiles"],
        canonical_smiles=result["canonical_smiles"],
        properties=PhysicochemicalProperties(**result["properties"]),
        pharmacokinetics=PharmacokineticProperties(**result["pharmacokinetics"]),
        toxicity=ToxicityProperties(**result["toxicity"]),
        warnings=result["warnings"]
    )
    
    return admet_response
