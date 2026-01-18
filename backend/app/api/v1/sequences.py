"""
Sequence analysis API endpoints.
"""
from fastapi import APIRouter, Depends, HTTPException, status, UploadFile, File
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy import select

from app.core.database import get_db
from app.core.security import get_current_active_user
from app.models.user import User
from app.models.sequence import Sequence
from app.services.sequence_service import SequenceService

router = APIRouter(prefix="/sequences", tags=["Sequences"])


@router.post("/import/fasta")
async def import_fasta(
    file: UploadFile = File(...),
    project_id: str = None,
    current_user: User = Depends(get_current_active_user),
    db: AsyncSession = Depends(get_db)
):
    """
    Import sequences from FASTA file.
    """
    content = await file.read()
    fasta_text = content.decode('utf-8')
    
    sequences = SequenceService.parse_fasta(fasta_text)
    
    if not sequences:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="No valid sequences found in FASTA file"
        )
    
    # Store sequences in database
    created_sequences = []
    for seq_data in sequences:
        is_valid, error = SequenceService.validate_protein_sequence(seq_data['sequence'])
        
        if not is_valid:
            continue  # Skip invalid sequences
        
        sequence = Sequence(
            project_id=project_id,
            source="upload",
            accession=seq_data['id'],
            sequence=seq_data['sequence'],
            metadata={
                'description': seq_data['description'],
                'length': seq_data['length']
            }
        )
        
        db.add(sequence)
        created_sequences.append(sequence)
    
    await db.commit()
    
    return {
        'imported': len(created_sequences),
        'total': len(sequences),
        'sequences': [{'id': str(s.id), 'accession': s.accession} for s in created_sequences]
    }


@router.post("/align/pairwise")
async def align_pairwise(
    seq1_id: str,
    seq2_id: str,
    current_user: User = Depends(get_current_active_user),
    db: AsyncSession = Depends(get_db)
):
    """
    Perform pairwise sequence alignment.
    """
    # Fetch sequences
    result1 = await db.execute(select(Sequence).where(Sequence.id == seq1_id))
    seq1 = result1.scalar_one_or_none()
    
    result2 = await db.execute(select(Sequence).where(Sequence.id == seq2_id))
    seq2 = result2.scalar_one_or_none()
    
    if not seq1 or not seq2:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="One or both sequences not found"
        )
    
    # Perform alignment
    alignment = SequenceService.pairwise_align_global(seq1.sequence, seq2.sequence)
    
    return alignment


@router.post("/properties")
async def calculate_properties(
    sequence_id: str,
    current_user: User = Depends(get_current_active_user),
    db: AsyncSession = Depends(get_db)
):
    """
    Calculate sequence properties.
    """
    result = await db.execute(select(Sequence).where(Sequence.id == sequence_id))
    sequence = result.scalar_one_or_none()
    
    if not sequence:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Sequence not found"
        )
    
    properties = SequenceService.calculate_sequence_properties(sequence.sequence)
    
    return properties
