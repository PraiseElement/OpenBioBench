"""
Database API endpoints for fetching from external biological databases.
"""
from fastapi import APIRouter, Depends, HTTPException, status, Query
from sqlalchemy.ext.asyncio import AsyncSession
from pydantic import BaseModel
from typing import List, Optional

from app.core.database import get_db
from app.core.security import get_current_active_user
from app.models.user import User
from app.services.external_db_service import ExternalDBService


router = APIRouter(prefix="/database", tags=["Database"])


# Response models
class UniProtSequence(BaseModel):
    accession: str
    header: str
    sequence: str
    length: int
    fasta: str
    metadata: dict


class UniProtSearchResult(BaseModel):
    accession: str
    entry_name: str
    protein_name: str
    gene_names: List[str]
    organism: str
    length: int


class PDBStructure(BaseModel):
    pdb_id: str
    pdb_content: str
    atom_count: int
    chains: List[str]
    metadata: dict


class PDBSearchResult(BaseModel):
    pdb_id: str
    title: str
    resolution: Optional[float]
    method: str


class AlphaFoldStructure(BaseModel):
    uniprot_id: str
    pdb_content: str
    model_version: str
    mean_plddt: Optional[float]
    created_date: str


class PubChemCompound(BaseModel):
    cid: int
    smiles: str
    iupac_name: str
    molecular_formula: str
    molecular_weight: float


# UniProt endpoints
@router.get("/uniprot/{accession}", response_model=UniProtSequence)
async def get_uniprot_sequence(
    accession: str,
    current_user: User = Depends(get_current_active_user)
):
    """
    Fetch protein sequence from UniProt by accession.
    
    Example accessions: P00533 (EGFR), P04637 (p53), Q9Y6K9
    """
    try:
        result = await ExternalDBService.fetch_uniprot_sequence(accession)
        return result
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to fetch from UniProt: {str(e)}")


@router.get("/uniprot/search/", response_model=List[UniProtSearchResult])
async def search_uniprot(
    query: str = Query(..., min_length=2),
    limit: int = Query(10, ge=1, le=50),
    current_user: User = Depends(get_current_active_user)
):
    """
    Search UniProt for proteins by name, gene, or organism.
    """
    try:
        results = await ExternalDBService.search_uniprot(query, limit)
        return results
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Search failed: {str(e)}")


# PDB endpoints
@router.get("/pdb/{pdb_id}", response_model=PDBStructure)
async def get_pdb_structure(
    pdb_id: str,
    current_user: User = Depends(get_current_active_user)
):
    """
    Fetch crystal structure from RCSB PDB by ID.
    
    Example IDs: 1HSG (HIV protease), 4HHB (hemoglobin), 6LU7 (SARS-CoV-2)
    """
    try:
        result = await ExternalDBService.fetch_pdb_structure(pdb_id)
        return result
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to fetch from PDB: {str(e)}")


@router.get("/pdb/search/", response_model=List[PDBSearchResult])
async def search_pdb(
    query: str = Query(..., min_length=2),
    limit: int = Query(10, ge=1, le=50),
    current_user: User = Depends(get_current_active_user)
):
    """
    Search RCSB PDB for structures.
    """
    try:
        results = await ExternalDBService.search_pdb(query, limit)
        return results
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Search failed: {str(e)}")


# AlphaFold endpoints
@router.get("/alphafold/{uniprot_id}", response_model=AlphaFoldStructure)
async def get_alphafold_structure(
    uniprot_id: str,
    current_user: User = Depends(get_current_active_user)
):
    """
    Fetch predicted structure from AlphaFold DB by UniProt accession.
    """
    try:
        result = await ExternalDBService.fetch_alphafold_structure(uniprot_id)
        return result
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to fetch from AlphaFold: {str(e)}")


# PubChem endpoints
@router.get("/pubchem/{cid}", response_model=PubChemCompound)
async def get_pubchem_compound(
    cid: int,
    current_user: User = Depends(get_current_active_user)
):
    """
    Fetch compound from PubChem by CID.
    
    Example CIDs: 2244 (Aspirin), 5090 (Diazepam), 679 (Dopamine)
    """
    try:
        result = await ExternalDBService.fetch_pubchem_compound(cid)
        return result
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to fetch from PubChem: {str(e)}")


@router.get("/pubchem/search/", response_model=List[PubChemCompound])
async def search_pubchem(
    query: str = Query(..., min_length=2),
    limit: int = Query(10, ge=1, le=20),
    current_user: User = Depends(get_current_active_user)
):
    """
    Search PubChem for compounds by name.
    """
    try:
        results = await ExternalDBService.search_pubchem(query, limit)
        return results
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Search failed: {str(e)}")
