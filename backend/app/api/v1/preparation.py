"""
Structure preparation API endpoints.
Handles protein and ligand preparation for docking.
"""
from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.ext.asyncio import AsyncSession
from pydantic import BaseModel
from typing import List, Optional, Dict

from app.core.database import get_db
from app.core.security import get_current_active_user
from app.models.user import User
from app.services.preparation_service import PreparationService


router = APIRouter(prefix="/preparation", tags=["Preparation"])


# Request/Response models
class CleanProteinRequest(BaseModel):
    pdb_content: str
    remove_water: bool = True
    remove_heteroatoms: bool = False
    keep_ligands: List[str] = []
    select_chains: Optional[List[str]] = None


class CleanProteinResponse(BaseModel):
    pdb_content: str
    stats: dict


class AddHydrogensRequest(BaseModel):
    pdb_content: str
    ph: float = 7.0


class AddHydrogensResponse(BaseModel):
    pdb_content: str
    hydrogens_present: int
    ph: float
    note: str


class PrepareLigandSmilesRequest(BaseModel):
    smiles: str
    optimize: bool = True
    num_conformers: int = 1


class PrepareLigandSdfRequest(BaseModel):
    sdf_content: str
    optimize: bool = True


class PrepareLigandResponse(BaseModel):
    smiles: str
    pdb_content: str
    mol_content: str
    num_atoms: int
    properties: dict
    num_conformers: Optional[int] = 1


class DockingBoxRequest(BaseModel):
    pdb_content: str
    ligand_residue: Optional[str] = None


class DockingBoxResponse(BaseModel):
    center: List[float]
    size: List[float]
    based_on: str


# Protein preparation endpoints
@router.post("/protein/clean", response_model=CleanProteinResponse)
async def clean_protein(
    request: CleanProteinRequest,
    current_user: User = Depends(get_current_active_user)
):
    """
    Clean a protein structure by removing waters and other unwanted atoms.
    """
    try:
        result = PreparationService.clean_protein(
            pdb_content=request.pdb_content,
            options={
                "remove_water": request.remove_water,
                "remove_heteroatoms": request.remove_heteroatoms,
                "keep_ligands": request.keep_ligands,
                "select_chains": request.select_chains
            }
        )
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Protein cleaning failed: {str(e)}")


@router.post("/protein/add_hydrogens", response_model=AddHydrogensResponse)
async def add_hydrogens(
    request: AddHydrogensRequest,
    current_user: User = Depends(get_current_active_user)
):
    """
    Add hydrogen atoms to a protein structure.
    """
    try:
        result = PreparationService.add_hydrogens_to_protein(
            pdb_content=request.pdb_content,
            ph=request.ph
        )
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Adding hydrogens failed: {str(e)}")


# Ligand preparation endpoints
@router.post("/ligand/from_smiles", response_model=PrepareLigandResponse)
async def prepare_ligand_smiles(
    request: PrepareLigandSmilesRequest,
    current_user: User = Depends(get_current_active_user)
):
    """
    Prepare a ligand from SMILES string for docking.
    Generates 3D coordinates and optimizes geometry.
    """
    try:
        result = PreparationService.prepare_ligand_from_smiles(
            smiles=request.smiles,
            options={
                "optimize": request.optimize,
                "num_conformers": request.num_conformers
            }
        )
        return result
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Ligand preparation failed: {str(e)}")


@router.post("/ligand/from_sdf", response_model=PrepareLigandResponse)
async def prepare_ligand_sdf(
    request: PrepareLigandSdfRequest,
    current_user: User = Depends(get_current_active_user)
):
    """
    Prepare a ligand from SDF file content.
    """
    try:
        result = PreparationService.prepare_ligand_from_sdf(
            sdf_content=request.sdf_content,
            options={"optimize": request.optimize}
        )
        return result
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Ligand preparation failed: {str(e)}")


# Utility endpoints
@router.post("/docking_box", response_model=DockingBoxResponse)
async def calculate_docking_box(
    request: DockingBoxRequest,
    current_user: User = Depends(get_current_active_user)
):
    """
    Calculate suggested docking box dimensions from a structure.
    """
    try:
        result = PreparationService.calculate_docking_box(
            pdb_content=request.pdb_content,
            ligand_residue=request.ligand_residue
        )
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Box calculation failed: {str(e)}")
