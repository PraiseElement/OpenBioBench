"""
Builder API endpoints for ligand design.
Handles 2D validation, 3D generation, and conformer search.
"""
import logging
import uuid
from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.ext.asyncio import AsyncSession
from rdkit import Chem

from app.core.database import get_db
from app.core.security import get_current_active_user
from app.models.user import User
from app.models.ligand import Ligand
from app.services.builder_service import BuilderService
from app.services.storage_service import storage
from app.schemas.builder import (
    SMILESValidateRequest,
    SMILESValidateResponse,
    Generate3DRequest,
    Generate3DResponse,
    ConformerRequest,
    ConformerResponse,
    ConformerInfo
)

logger = logging.getLogger(__name__)
router = APIRouter(prefix="/builder", tags=["Ligand Builder"])


@router.post("/validate", response_model=SMILESValidateResponse)
async def validate_smiles(
    request: SMILESValidateRequest,
    current_user: User = Depends(get_current_active_user)
):
    """
    Validate a SMILES string and return canonical form with basic properties.
    """
    is_valid, canonical, error = BuilderService.validate_smiles(request.smiles)
    
    if not is_valid:
        return SMILESValidateResponse(
            valid=False,
            error=error
        )
    
    # Calculate basic properties
    mol = Chem.MolFromSmiles(canonical)
    props = BuilderService.calculate_properties(mol)
    
    return SMILESValidateResponse(
        valid=True,
        canonical_smiles=canonical,
        molecular_formula=props["molecular_formula"],
        molecular_weight=props["molecular_weight"],
        error=None
    )


@router.post("/2d/to_3d", response_model=Generate3DResponse)
async def generate_3d(
    request: Generate3DRequest,
    current_user: User = Depends(get_current_active_user),
    db: AsyncSession = Depends(get_db)
):
    """
    Generate 3D conformers from a SMILES string.
    
    Creates a ligand entry in the database and stores conformer files.
    """
    # Validate SMILES
    is_valid, canonical, error = BuilderService.validate_smiles(request.smiles)
    if not is_valid:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Invalid SMILES: {error}"
        )
    
    # Generate conformers
    conformers = BuilderService.generate_conformers(
        canonical,
        num_conformers=request.num_conformers,
        method=request.method,
        force_field=request.force_field if request.minimize else None,
        random_seed=42
    )
    
    if not conformers:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to generate 3D conformers"
        )
    
    # Calculate properties
    mol = Chem.MolFromSmiles(canonical)
    properties = BuilderService.calculate_properties(mol)
    inchi, inchi_key = BuilderService.calculate_inchi(mol)
    
    # Create ligand entry
    ligand = Ligand(
        project_id=request.project_id,
        source="builder",
        smiles=canonical,
        inchi=inchi,
        inchi_key=inchi_key,
        name=request.name,
        properties=properties
    )
    
    db.add(ligand)
    await db.commit()
    await db.refresh(ligand)
    
    # Store conformer files
    conformer_info_list = []
    lowest_energy_id = 0
    lowest_energy = float('inf')
    
    for idx, (mol_conf, energy, conf_id) in enumerate(conformers):
        # Write MOL2 file
        mol2_content = Chem.MolToMolBlock(mol_conf, confId=conf_id)
        file_url = storage.upload_text(
            mol2_content,
            f"ligand_{ligand.id}_conf_{idx}.mol2",
            folder="ligands"
        )
        
        # Calculate RMSD to first conformer
        rmsd = None
        if idx > 0:
            from rdkit.Chem import AllChem
            rmsd = AllChem.GetConformerRMS(mol_conf, conf_id, conformers[0][2])
        
        conformer_info_list.append(
            ConformerInfo(
                conformer_id=idx,
                energy=round(energy, 2),
                rmsd=round(rmsd, 3) if rmsd else None,
                file_url=file_url
            )
        )
        
        if energy < lowest_energy:
            lowest_energy = energy
            lowest_energy_id = idx
    
    # Update ligand with lowest energy conformer
    ligand.mol_file_url = conformer_info_list[lowest_energy_id].file_url
    await db.commit()
    
    return Generate3DResponse(
        ligand_id=ligand.id,
        smiles=request.smiles,
        canonical_smiles=canonical,
        conformers=conformer_info_list,
        lowest_energy_conformer_id=lowest_energy_id,
        properties=properties
    )


@router.post("/3d/optimize")
async def optimize_3d(
    ligand_id: uuid.UUID,
    force_field: str = "mmff94",
    current_user: User = Depends(get_current_active_user),
    db: AsyncSession = Depends(get_db)
):
    """
    Energy minimize an existing 3D ligand structure.
    """
    # Implementation would load ligand, minimize, and update
    raise HTTPException(
        status_code=status.HTTP_501_NOT_IMPLEMENTED,
        detail="Endpoint under development"
    )
