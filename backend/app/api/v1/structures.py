"""
Structure analysis API endpoints.
"""
from fastapi import APIRouter, Depends, HTTPException, status, UploadFile, File
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy import select

from app.core.database import get_db
from app.core.security import get_current_active_user
from app.models.user import User
from app.models.structure import Structure
from app.services.structure_service import StructureService
from app.services.storage_service import storage

router = APIRouter(prefix="/structures", tags=["Structures"])


@router.post("/import/pdb")
async def import_pdb(
    file: UploadFile = File(...),
    project_id: str = None,
    pdb_id: str = None,
    current_user: User = Depends(get_current_active_user),
    db: AsyncSession = Depends(get_db)
):
    """
    Import protein structure from PDB file.
    """
    content = await file.read()
    pdb_text = content.decode('utf-8')
    
    # Parse structure
    info = StructureService.parse_pdb(pdb_text)
    
    if not info:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Failed to parse PDB file"
        )
    
    # Validate structure
    is_valid, issues = StructureService.validate_structure(pdb_text)
    
    # Upload to storage
    file_url = storage.upload_text(pdb_text, file.filename or 'structure.pdb', 'structures')
    
    # Create structure entry
    structure = Structure(
        project_id=project_id,
        source="upload",
        pdb_id=pdb_id,
        file_url=file_url,
        file_format="pdb",
        method="upload",
        chains=[c['chain_id'] for c in info['chains']],
        num_residues=info['num_residues'],
        metadata={
            'validation_issues': issues,
            'is_valid': is_valid,
            'sequence': info.get('sequence', '')
        }
    )
    
    db.add(structure)
    await db.commit()
    await db.refresh(structure)
    
    return {
        'structure_id': str(structure.id),
        'chains': info['chains'],
        'num_residues': info['num_residues'],
        'is_valid': is_valid,
        'issues': issues
    }


@router.post("/{structure_id}/pockets")
async def detect_pockets(
    structure_id: str,
    current_user: User = Depends(get_current_active_user),
    db: AsyncSession = Depends(get_db)
):
    """
    Detect binding pockets in a structure.
    """
    result = await db.execute(select(Structure).where(Structure.id == structure_id))
    structure = result.scalar_one_or_none()
    
    if not structure:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Structure not found"
        )
    
    # Download structure file
    pdb_content = storage.download_file(structure.file_url)
    if not pdb_content:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Failed to load structure file"
        )
    
    # Detect pockets
    pockets = StructureService.detect_binding_pockets(pdb_content.decode('utf-8'))
    
    return {
        'structure_id': str(structure.id),
        'pockets': pockets,
        'num_pockets': len(pockets)
    }


@router.get("/{structure_id}")
async def get_structure(
    structure_id: str,
    current_user: User = Depends(get_current_active_user),
    db: AsyncSession = Depends(get_db)
):
    """
    Get structure information.
    """
    result = await db.execute(select(Structure).where(Structure.id == structure_id))
    structure = result.scalar_one_or_none()
    
    if not structure:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Structure not found"
        )
    
    return structure
