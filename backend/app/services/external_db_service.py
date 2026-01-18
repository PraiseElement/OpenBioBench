"""
External Database Service.
Provides unified access to biological databases: UniProt, PDB, AlphaFold, PubChem.
"""
import httpx
import logging
from typing import Dict, Optional, List
from io import StringIO

logger = logging.getLogger(__name__)

# API Base URLs
UNIPROT_API = "https://rest.uniprot.org"
PDB_API = "https://data.rcsb.org/rest/v1"
PDB_FILES = "https://files.rcsb.org/download"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api"
PUBCHEM_API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
EBI_TOOLS = "https://www.ebi.ac.uk/Tools/services/rest"


class ExternalDBService:
    """Service for fetching data from external biological databases."""
    
    @staticmethod
    async def fetch_uniprot_sequence(accession: str) -> Dict:
        """
        Fetch protein sequence and metadata from UniProt.
        
        Args:
            accession: UniProt accession (e.g., P00533, Q9Y6K9)
            
        Returns:
            Dictionary with sequence and metadata
        """
        async with httpx.AsyncClient(timeout=30.0) as client:
            # Fetch FASTA format
            fasta_url = f"{UNIPROT_API}/uniprotkb/{accession}.fasta"
            fasta_response = await client.get(fasta_url)
            
            if fasta_response.status_code == 404:
                raise ValueError(f"UniProt accession {accession} not found")
            fasta_response.raise_for_status()
            
            fasta_content = fasta_response.text
            
            # Parse FASTA
            lines = fasta_content.strip().split('\n')
            header = lines[0][1:] if lines[0].startswith('>') else ""
            sequence = ''.join(lines[1:])
            
            # Fetch JSON metadata
            json_url = f"{UNIPROT_API}/uniprotkb/{accession}"
            json_response = await client.get(json_url, headers={"Accept": "application/json"})
            
            metadata = {}
            if json_response.status_code == 200:
                try:
                    data = json_response.json()
                    metadata = {
                        "organism": data.get("organism", {}).get("scientificName", "Unknown"),
                        "protein_name": data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", ""),
                        "gene_names": [g.get("geneName", {}).get("value", "") for g in data.get("genes", [])],
                        "function": "",
                    }
                    # Extract function from comments
                    for comment in data.get("comments", []):
                        if comment.get("commentType") == "FUNCTION":
                            texts = comment.get("texts", [])
                            if texts:
                                metadata["function"] = texts[0].get("value", "")
                except Exception as e:
                    logger.warning(f"Failed to parse UniProt JSON: {e}")
            
            return {
                "accession": accession,
                "header": header,
                "sequence": sequence,
                "length": len(sequence),
                "fasta": fasta_content,
                "metadata": metadata
            }
    
    @staticmethod
    async def search_uniprot(query: str, limit: int = 10) -> List[Dict]:
        """
        Search UniProt for proteins matching query.
        
        Args:
            query: Search query (gene name, protein name, organism, etc.)
            limit: Maximum results to return
            
        Returns:
            List of matching entries with basic info
        """
        async with httpx.AsyncClient(timeout=30.0) as client:
            url = f"{UNIPROT_API}/uniprotkb/search"
            params = {
                "query": query,
                "format": "json",
                "size": limit,
                "fields": "accession,id,protein_name,gene_names,organism_name,length"
            }
            
            response = await client.get(url, params=params)
            response.raise_for_status()
            
            data = response.json()
            results = []
            
            for entry in data.get("results", []):
                results.append({
                    "accession": entry.get("primaryAccession", ""),
                    "entry_name": entry.get("uniProtkbId", ""),
                    "protein_name": entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", ""),
                    "gene_names": [g.get("geneName", {}).get("value", "") for g in entry.get("genes", [])],
                    "organism": entry.get("organism", {}).get("scientificName", ""),
                    "length": entry.get("sequence", {}).get("length", 0)
                })
            
            return results
    
    @staticmethod
    async def fetch_pdb_structure(pdb_id: str) -> Dict:
        """
        Fetch PDB structure file and metadata from RCSB PDB.
        
        Args:
            pdb_id: PDB ID (e.g., 1HSG, 4HHB)
            
        Returns:
            Dictionary with PDB content and metadata
        """
        pdb_id = pdb_id.upper().strip()
        
        async with httpx.AsyncClient(timeout=60.0) as client:
            # Fetch PDB file
            pdb_url = f"{PDB_FILES}/{pdb_id}.pdb"
            pdb_response = await client.get(pdb_url)
            
            if pdb_response.status_code == 404:
                raise ValueError(f"PDB ID {pdb_id} not found")
            pdb_response.raise_for_status()
            
            pdb_content = pdb_response.text
            
            # Fetch metadata
            meta_url = f"{PDB_API}/core/entry/{pdb_id}"
            meta_response = await client.get(meta_url)
            
            metadata = {}
            if meta_response.status_code == 200:
                try:
                    data = meta_response.json()
                    metadata = {
                        "title": data.get("struct", {}).get("title", ""),
                        "resolution": data.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0],
                        "method": data.get("exptl", [{}])[0].get("method", ""),
                        "release_date": data.get("rcsb_accession_info", {}).get("initial_release_date", ""),
                        "organism": "",
                    }
                except Exception as e:
                    logger.warning(f"Failed to parse PDB metadata: {e}")
            
            # Count atoms and chains
            atom_count = pdb_content.count("\nATOM")
            chains = set()
            for line in pdb_content.split('\n'):
                if line.startswith('ATOM') and len(line) > 21:
                    chains.add(line[21])
            
            return {
                "pdb_id": pdb_id,
                "pdb_content": pdb_content,
                "atom_count": atom_count,
                "chains": list(chains),
                "metadata": metadata
            }
    
    @staticmethod
    async def search_pdb(query: str, limit: int = 10) -> List[Dict]:
        """
        Search RCSB PDB for structures.
        
        Args:
            query: Search query
            limit: Maximum results
            
        Returns:
            List of matching structures
        """
        async with httpx.AsyncClient(timeout=30.0) as client:
            url = "https://search.rcsb.org/rcsbsearch/v2/query"
            
            search_request = {
                "query": {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {
                        "value": query
                    }
                },
                "return_type": "entry",
                "request_options": {
                    "paginate": {
                        "start": 0,
                        "rows": limit
                    }
                }
            }
            
            response = await client.post(url, json=search_request)
            
            if response.status_code != 200:
                return []
            
            data = response.json()
            results = []
            
            for hit in data.get("result_set", []):
                pdb_id = hit.get("identifier", "")
                # Fetch basic info for each hit
                meta_url = f"{PDB_API}/core/entry/{pdb_id}"
                meta_response = await client.get(meta_url)
                
                if meta_response.status_code == 200:
                    meta = meta_response.json()
                    results.append({
                        "pdb_id": pdb_id,
                        "title": meta.get("struct", {}).get("title", ""),
                        "resolution": meta.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0],
                        "method": meta.get("exptl", [{}])[0].get("method", "")
                    })
            
            return results
    
    @staticmethod
    async def fetch_alphafold_structure(uniprot_id: str) -> Dict:
        """
        Fetch AlphaFold predicted structure for a UniProt accession.
        
        Args:
            uniprot_id: UniProt accession
            
        Returns:
            Dictionary with PDB content and pLDDT scores
        """
        async with httpx.AsyncClient(timeout=60.0) as client:
            # Get AlphaFold entry info
            url = f"{ALPHAFOLD_API}/prediction/{uniprot_id}"
            response = await client.get(url)
            
            if response.status_code == 404:
                raise ValueError(f"AlphaFold structure not available for {uniprot_id}")
            response.raise_for_status()
            
            data = response.json()
            if not data:
                raise ValueError(f"No AlphaFold data for {uniprot_id}")
            
            entry = data[0] if isinstance(data, list) else data
            
            # Download PDB file
            pdb_url = entry.get("pdbUrl", "")
            if not pdb_url:
                raise ValueError("No PDB URL in AlphaFold response")
            
            pdb_response = await client.get(pdb_url)
            pdb_response.raise_for_status()
            
            return {
                "uniprot_id": uniprot_id,
                "pdb_content": pdb_response.text,
                "model_version": entry.get("latestVersion", "unknown"),
                "mean_plddt": entry.get("meanPlddt", None),
                "created_date": entry.get("modelCreatedDate", "")
            }
    
    @staticmethod
    async def fetch_pubchem_compound(cid: int) -> Dict:
        """
        Fetch compound data from PubChem.
        
        Args:
            cid: PubChem Compound ID
            
        Returns:
            Dictionary with compound info and SMILES
        """
        async with httpx.AsyncClient(timeout=30.0) as client:
            # Get compound properties
            url = f"{PUBCHEM_API}/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IUPACName/JSON"
            response = await client.get(url)
            
            if response.status_code == 404:
                raise ValueError(f"PubChem CID {cid} not found")
            response.raise_for_status()
            
            data = response.json()
            props = data.get("PropertyTable", {}).get("Properties", [{}])[0]
            
            return {
                "cid": cid,
                "smiles": props.get("CanonicalSMILES", ""),
                "iupac_name": props.get("IUPACName", ""),
                "molecular_formula": props.get("MolecularFormula", ""),
                "molecular_weight": props.get("MolecularWeight", 0)
            }
    
    @staticmethod
    async def search_pubchem(query: str, limit: int = 10) -> List[Dict]:
        """
        Search PubChem for compounds by name.
        
        Args:
            query: Compound name or partial name
            limit: Maximum results
            
        Returns:
            List of matching compounds
        """
        async with httpx.AsyncClient(timeout=30.0) as client:
            # Search by name
            url = f"{PUBCHEM_API}/compound/name/{query}/cids/JSON"
            response = await client.get(url)
            
            if response.status_code != 200:
                return []
            
            data = response.json()
            cids = data.get("IdentifierList", {}).get("CID", [])[:limit]
            
            results = []
            for cid in cids:
                try:
                    compound = await ExternalDBService.fetch_pubchem_compound(cid)
                    results.append(compound)
                except Exception:
                    continue
            
            return results
