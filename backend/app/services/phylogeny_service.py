"""
Phylogenetic analysis service.
Builds phylogenetic trees from sequence alignments.
"""
from typing import List, Dict, Tuple
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO


class PhylogenyService:
    """Service for phylogenetic tree construction."""

    @staticmethod
    def parse_aligned_fasta(fasta_content: str) -> List[SeqRecord]:
        """Parse aligned FASTA sequences."""
        records = []
        current_id = None
        current_seq = []

        for line in fasta_content.strip().split('\n'):
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    seq = ''.join(current_seq).upper()
                    records.append(SeqRecord(Seq(seq), id=current_id, description=""))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_id is not None:
            seq = ''.join(current_seq).upper()
            records.append(SeqRecord(Seq(seq), id=current_id, description=""))

        return records

    @staticmethod
    def build_tree(alignment_fasta: str, method: str = "nj") -> Dict:
        """
        Build phylogenetic tree from sequence alignment.
        
        Args:
            alignment_fasta: Aligned sequences in FASTA format
            method: Tree construction method ('nj' or 'upgma')
            
        Returns:
            Dictionary with tree data
        """
        # Parse sequences
        records = PhylogenyService.parse_aligned_fasta(alignment_fasta)
        
        if len(records) < 3:
            raise ValueError("At least 3 sequences required for phylogenetic analysis")
        
        if len(records) > 500:
            raise ValueError("Maximum 500 sequences supported")
        
        # Pad sequences to same length
        max_len = max(len(str(r.seq)) for r in records)
        for r in records:
            if len(str(r.seq)) < max_len:
                r.seq = Seq(str(r.seq) + '-' * (max_len - len(str(r.seq))))
        
        # Create multiple sequence alignment
        alignment = MultipleSeqAlignment(records)
        
        # Calculate distance matrix
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        
        # Build tree
        constructor = DistanceTreeConstructor()
        
        if method.lower() == 'upgma':
            tree = constructor.upgma(dm)
        else:  # Default to NJ
            tree = constructor.nj(dm)
        
        # Convert tree to Newick format
        output = StringIO()
        Phylo.write(tree, output, 'newick')
        newick = output.getvalue().strip()
        
        # Generate ASCII tree representation
        ascii_output = StringIO()
        Phylo.draw_ascii(tree, ascii_output)
        ascii_tree = ascii_output.getvalue()
        
        # Calculate total branch length
        total_length = sum(
            clade.branch_length or 0 
            for clade in tree.find_clades()
        )
        
        return {
            'newick': newick,
            'ascii_tree': ascii_tree,
            'num_taxa': len(records),
            'method': method.upper(),
            'total_branch_length': round(total_length, 4)
        }

    @staticmethod
    def calculate_distance_matrix(alignment_fasta: str) -> Dict:
        """Calculate pairwise distance matrix from alignment."""
        records = PhylogenyService.parse_aligned_fasta(alignment_fasta)
        
        # Pad sequences
        max_len = max(len(str(r.seq)) for r in records)
        for r in records:
            if len(str(r.seq)) < max_len:
                r.seq = Seq(str(r.seq) + '-' * (max_len - len(str(r.seq))))
        
        alignment = MultipleSeqAlignment(records)
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        
        # Convert to dict format
        names = [r.id for r in records]
        matrix = []
        for i, row in enumerate(dm.matrix):
            matrix.append(row)
        
        return {
            'names': names,
            'matrix': matrix
        }
