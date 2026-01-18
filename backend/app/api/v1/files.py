"""
File upload API endpoints.
Handles PDB, PDBQT, SDF, MOL2 file uploads for molecular analysis.
"""
from fastapi import APIRouter, Depends, HTTPException, status, UploadFile, File, Form
from sqlalchemy.ext.asyncio import AsyncSession
from pydantic import BaseModel
from typing import Optional, List
import tempfile
import os

from app.core.database import get_db
from app.core.security import get_current_active_user
from app.models.user import User


router = APIRouter(prefix="/files", tags=["Files"])


class FileUploadResponse(BaseModel):
    filename: str
    file_type: str
    size_bytes: int
    content: str
    validation: dict


def validate_pdb_content(content: str) -> dict:
    """Validate PDB file content."""
    lines = content.split('\n')
    atom_count = sum(1 for line in lines if line.startswith('ATOM'))
    hetatm_count = sum(1 for line in lines if line.startswith('HETATM'))
    
    chains = set()
    for line in lines:
        if (line.startswith('ATOM') or line.startswith('HETATM')) and len(line) > 21:
            chains.add(line[21])
    
    return {
        "valid": atom_count > 0 or hetatm_count > 0,
        "atom_count": atom_count,
        "hetatm_count": hetatm_count,
        "chains": list(chains),
        "line_count": len(lines)
    }


def validate_sdf_content(content: str) -> dict:
    """Validate SDF file content."""
    molecules = content.split('$$$$')
    mol_count = len([m for m in molecules if m.strip()])
    
    return {
        "valid": mol_count > 0,
        "molecule_count": mol_count,
        "size_bytes": len(content)
    }


def validate_mol2_content(content: str) -> dict:
    """Validate MOL2 file content."""
    has_molecule = '@<TRIPOS>MOLECULE' in content
    has_atoms = '@<TRIPOS>ATOM' in content
    
    atom_count = 0
    if has_atoms:
        in_atom_block = False
        for line in content.split('\n'):
            if '@<TRIPOS>ATOM' in line:
                in_atom_block = True
                continue
            if '@<TRIPOS>' in line and in_atom_block:
                break
            if in_atom_block and line.strip():
                atom_count += 1
    
    return {
        "valid": has_molecule and has_atoms,
        "atom_count": atom_count,
        "has_molecule_record": has_molecule,
        "has_atom_block": has_atoms
    }


@router.post("/upload/pdb", response_model=FileUploadResponse)
async def upload_pdb_file(
    file: UploadFile = File(...),
    current_user: User = Depends(get_current_active_user)
):
    """
    Upload a PDB file for molecular analysis.
    """
    if not file.filename.lower().endswith('.pdb'):
        raise HTTPException(
            status_code=400,
            detail="File must have .pdb extension"
        )
    
    content = await file.read()
    try:
        content_str = content.decode('utf-8')
    except UnicodeDecodeError:
        raise HTTPException(status_code=400, detail="Invalid file encoding")
    
    validation = validate_pdb_content(content_str)
    
    if not validation["valid"]:
        raise HTTPException(
            status_code=400,
            detail="Invalid PDB file: no ATOM records found"
        )
    
    return FileUploadResponse(
        filename=file.filename,
        file_type="pdb",
        size_bytes=len(content),
        content=content_str,
        validation=validation
    )


@router.post("/upload/pdbqt", response_model=FileUploadResponse)
async def upload_pdbqt_file(
    file: UploadFile = File(...),
    current_user: User = Depends(get_current_active_user)
):
    """
    Upload a PDBQT file (AutoDock format) for docking.
    """
    if not file.filename.lower().endswith('.pdbqt'):
        raise HTTPException(
            status_code=400,
            detail="File must have .pdbqt extension"
        )
    
    content = await file.read()
    try:
        content_str = content.decode('utf-8')
    except UnicodeDecodeError:
        raise HTTPException(status_code=400, detail="Invalid file encoding")
    
    # PDBQT validation (similar to PDB but with charge/type columns)
    validation = validate_pdb_content(content_str)
    validation["file_type"] = "pdbqt"
    
    return FileUploadResponse(
        filename=file.filename,
        file_type="pdbqt",
        size_bytes=len(content),
        content=content_str,
        validation=validation
    )


@router.post("/upload/sdf", response_model=FileUploadResponse)
async def upload_sdf_file(
    file: UploadFile = File(...),
    current_user: User = Depends(get_current_active_user)
):
    """
    Upload an SDF file containing one or more molecules.
    """
    if not file.filename.lower().endswith('.sdf'):
        raise HTTPException(
            status_code=400,
            detail="File must have .sdf extension"
        )
    
    content = await file.read()
    try:
        content_str = content.decode('utf-8')
    except UnicodeDecodeError:
        raise HTTPException(status_code=400, detail="Invalid file encoding")
    
    validation = validate_sdf_content(content_str)
    
    if not validation["valid"]:
        raise HTTPException(
            status_code=400,
            detail="Invalid SDF file: no molecules found"
        )
    
    return FileUploadResponse(
        filename=file.filename,
        file_type="sdf",
        size_bytes=len(content),
        content=content_str,
        validation=validation
    )


@router.post("/upload/mol2", response_model=FileUploadResponse)
async def upload_mol2_file(
    file: UploadFile = File(...),
    current_user: User = Depends(get_current_active_user)
):
    """
    Upload a MOL2 file for molecular analysis.
    """
    if not file.filename.lower().endswith('.mol2'):
        raise HTTPException(
            status_code=400,
            detail="File must have .mol2 extension"
        )
    
    content = await file.read()
    try:
        content_str = content.decode('utf-8')
    except UnicodeDecodeError:
        raise HTTPException(status_code=400, detail="Invalid file encoding")
    
    validation = validate_mol2_content(content_str)
    
    if not validation["valid"]:
        raise HTTPException(
            status_code=400,
            detail="Invalid MOL2 file: missing required records"
        )
    
    return FileUploadResponse(
        filename=file.filename,
        file_type="mol2",
        size_bytes=len(content),
        content=content_str,
        validation=validation
    )
