"""
Sequence Analysis Service.
Handles sequence operations, alignment, and analysis.
"""
import logging
from typing import List, Dict, Optional, Tuple
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import tempfile
import subprocess

logger = logging.getLogger(__name__)


class SequenceService:
    """Service for biological sequence analysis."""
    
    @staticmethod
    def parse_fasta(fasta_content: str) -> List[Dict]:
        """
        Parse FASTA format sequences.
        
        Args:
            fasta_content: FASTA formatted text
            
        Returns:
            List of sequence dictionaries
        """
        sequences = []
        
        try:
            for record in SeqIO.parse(StringIO(fasta_content), "fasta"):
                sequences.append({
                    'id': record.id,
                    'description': record.description,
                    'sequence': str(record.seq),
                    'length': len(record.seq)
                })
            
            return sequences
            
        except Exception as e:
            logger.error(f"FASTA parsing failed: {e}")
            return []
    
    @staticmethod
    def validate_protein_sequence(sequence: str) -> Tuple[bool, str]:
        """
        Validate protein sequence.
        
        Args:
            sequence: Amino acid sequence
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        # Standard amino acid codes
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        ambiguous = set('BZX')  # Allowed ambiguity codes
        
        sequence = sequence.upper().replace(' ', '').replace('\n', '')
        
        if len(sequence) == 0:
            return False, "Empty sequence"
        
        if len(sequence) < 20:
            return True, "Warning: Very short sequence (< 20 residues)"
        
        if len(sequence) > 50000:
            return False, "Sequence too long (> 50,000 residues)"
        
        # Check for invalid characters
        invalid_chars = set(sequence) - valid_aa - ambiguous
        if invalid_chars:
            return False, f"Invalid amino acids: {', '.join(invalid_chars)}"
        
        # Check for high ambiguity
        ambig_count = sum(1 for aa in sequence if aa in ambiguous)
        if ambig_count / len(sequence) > 0.1:
            return True, f"Warning: High ambiguity ({ambig_count / len(sequence) * 100:.1f}%)"
        
        return True, ""
    
    @staticmethod
    def pairwise_align_global(seq1: str, seq2: str) -> Dict:
        """
        Perform global pairwise alignment (Needleman-Wunsch).
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            
        Returns:
            Alignment result dictionary
        """
        try:
            aligner = Align.PairwiseAligner()
            aligner.mode = 'global'
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -2
            aligner.extend_gap_score = -0.5
            
            alignments = aligner.align(seq1, seq2)
            best_alignment = alignments[0]
            
            # Calculate statistics
            aligned_seq1 = str(best_alignment[0])
            aligned_seq2 = str(best_alignment[1])
            
            matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b and a != '-')
            gaps = aligned_seq1.count('-') + aligned_seq2.count('-')
            identity = matches / max(len(seq1), len(seq2)) * 100
            
            return {
                'aligned_seq1': aligned_seq1,
                'aligned_seq2': aligned_seq2,
                'score': best_alignment.score,
                'identity_percent': identity,
                'matches': matches,
                'gaps': gaps,
                'length': len(aligned_seq1)
            }
            
        except Exception as e:
            logger.error(f"Pairwise alignment failed: {e}")
            return {}
    
    @staticmethod
    def calculate_sequence_properties(sequence: str) -> Dict:
        """
        Calculate sequence properties.
        
        Args:
            sequence: Amino acid sequence
            
        Returns:
            Dictionary of properties
        """
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
        
        try:
            protein = ProteinAnalysis(sequence)
            
            properties = {
                'length': len(sequence),
                'molecular_weight': protein.molecular_weight(),
                'aromaticity': protein.aromaticity(),
                'instability_index': protein.instability_index(),
                'isoelectric_point': protein.isoelectric_point(),
                'gravy': protein.gravy(),  # Hydrophobicity
                'amino_acid_composition': protein.get_amino_acids_percent()
            }
            
            return properties
            
        except Exception as e:
            logger.error(f"Property calculation failed: {e}")
            return {}
    
    @staticmethod
    def run_mafft_alignment(sequences: List[Tuple[str, str]], method: str = "auto") -> str:
        """
        Run MAFFT multiple sequence alignment.
        
        Args:
            sequences: List of (id, sequence) tuples
            method: MAFFT method (auto, linsi, einsi, ginsi)
            
        Returns:
            Aligned sequences in FASTA format
        """
        # In production, this runs MAFFT in Docker container
        # Placeholder implementation
        
        logger.warning("MAFFT execution is placeholder - Docker container not available")
        
        # Return input sequences as "aligned" (placeholder)
        fasta_output = ""
        for seq_id, sequence in sequences:
            fasta_output += f">{seq_id}\n{sequence}\n"
        
        return fasta_output
    
    @staticmethod
    def calculate_conservation_scores(alignment: List[str]) -> List[float]:
        """
        Calculate Shannon entropy-based conservation scores.
        
        Args:
            alignment: List of aligned sequences
            
        Returns:
            List of conservation scores (0-1, higher = more conserved)
        """
        import math
        from collections import Counter
        
        if not alignment:
            return []
        
        alignment_length = len(alignment[0])
        scores = []
        
        for pos in range(alignment_length):
            # Get residues at this position
            residues = [seq[pos] for seq in alignment if pos < len(seq)]
            
            # Count occurrences
            counts = Counter(residues)
            total = len(residues)
            
            # Calculate Shannon entropy
            entropy = 0
            for count in counts.values():
                if count > 0:
                    p = count / total
                    entropy -= p * math.log2(p)
            
            # Normalize to 0-1 (max entropy for 20 amino acids is log2(20))
            max_entropy = math.log2(20)
            conservation = 1 - (entropy / max_entropy) if max_entropy > 0 else 0
            
            scores.append(conservation)
        
        return scores
