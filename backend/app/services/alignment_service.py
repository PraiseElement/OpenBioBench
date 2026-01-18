"""
Sequence alignment service.
Provides multiple sequence alignment with conservation analysis.
"""
from typing import List, Dict, Tuple
from Bio import pairwise2
from Bio.Seq import Seq
from collections import Counter
import math


class AlignmentService:
    """Service for sequence alignment operations."""

    @staticmethod
    def parse_fasta(fasta_content: str) -> List[Dict]:
        """Parse FASTA format content into list of sequences."""
        sequences = []
        current_id = None
        current_seq = []

        for line in fasta_content.strip().split('\n'):
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences.append({
                        'id': current_id,
                        'sequence': ''.join(current_seq)
                    })
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())

        if current_id is not None:
            sequences.append({
                'id': current_id,
                'sequence': ''.join(current_seq)
            })

        return sequences

    @staticmethod
    def perform_multiple_alignment(fasta_content: str) -> Dict:
        """
        Perform multiple sequence alignment.
        Uses progressive alignment approach with pairwise alignments.
        """
        # Parse sequences
        sequences = AlignmentService.parse_fasta(fasta_content)
        
        if len(sequences) < 2:
            raise ValueError("At least 2 sequences required for alignment")
        
        if len(sequences) > 100:
            raise ValueError("Maximum 100 sequences supported")

        # For MVP, use a simplified progressive alignment
        # Start with first sequence as reference
        aligned = []
        ref_seq = sequences[0]['sequence']
        
        # Align each sequence to reference
        for seq_data in sequences:
            if seq_data['id'] == sequences[0]['id']:
                aligned.append({
                    'id': seq_data['id'],
                    'sequence': ref_seq
                })
            else:
                # Perform pairwise alignment
                alignments = pairwise2.align.globalxx(
                    ref_seq, 
                    seq_data['sequence'],
                    one_alignment_only=True
                )
                if alignments:
                    aln = alignments[0]
                    aligned.append({
                        'id': seq_data['id'],
                        'sequence': aln.seqB
                    })
                    # Update reference to aligned version
                    ref_seq = aln.seqA

        # Pad sequences to same length
        max_len = max(len(s['sequence']) for s in aligned)
        for s in aligned:
            if len(s['sequence']) < max_len:
                s['sequence'] = s['sequence'] + '-' * (max_len - len(s['sequence']))

        # Calculate statistics
        alignment_length = max_len
        
        # Calculate conservation per position
        conservation = []
        identity_count = 0
        gap_count = 0
        
        for i in range(alignment_length):
            column = [s['sequence'][i] if i < len(s['sequence']) else '-' for s in aligned]
            
            # Count gaps
            gaps = column.count('-')
            gap_count += gaps
            
            # Calculate conservation (Shannon entropy-based)
            non_gap = [c for c in column if c != '-']
            if non_gap:
                counter = Counter(non_gap)
                total = len(non_gap)
                
                # Check if all same
                if len(counter) == 1:
                    conservation.append(1.0)
                    identity_count += 1
                else:
                    # Shannon entropy
                    entropy = 0
                    for count in counter.values():
                        p = count / total
                        entropy -= p * math.log2(p)
                    # Normalize: max entropy for amino acids is log2(20) ≈ 4.32
                    max_entropy = math.log2(min(20, total))
                    cons_score = 1 - (entropy / max_entropy) if max_entropy > 0 else 1
                    conservation.append(max(0, min(1, cons_score)))
            else:
                conservation.append(0.0)

        total_positions = alignment_length * len(aligned)
        
        return {
            'aligned_sequences': aligned,
            'num_sequences': len(aligned),
            'alignment_length': alignment_length,
            'percent_identity': identity_count / alignment_length if alignment_length > 0 else 0,
            'gap_percentage': gap_count / total_positions if total_positions > 0 else 0,
            'conservation': conservation
        }

    @staticmethod
    def pairwise_align_global(seq1: str, seq2: str) -> Dict:
        """Perform global pairwise alignment."""
        alignments = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)
        
        if not alignments:
            raise ValueError("Alignment failed")
        
        aln = alignments[0]
        
        # Calculate identity
        matches = sum(1 for a, b in zip(aln.seqA, aln.seqB) if a == b and a != '-')
        length = len(aln.seqA)
        
        return {
            'aligned_seq1': aln.seqA,
            'aligned_seq2': aln.seqB,
            'score': aln.score,
            'identity': matches / length if length > 0 else 0,
            'length': length
        }
