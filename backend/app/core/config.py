"""
Core configuration management using Pydantic Settings.
Loads environment variables and provides typed configuration access.
"""
from typing import List
from pydantic_settings import BaseSettings, SettingsConfigDict
from pydantic import AnyHttpUrl, PostgresDsn


class Settings(BaseSettings):
    """Application settings loaded from environment variables."""
    
    model_config = SettingsConfigDict(
        env_file=".env",
        case_sensitive=True,
        extra="ignore"
    )
    
    # Application
    APP_NAME: str = "OpenBioBench"
    APP_VERSION: str = "1.0.0"
    DEBUG: bool = False
    SECRET_KEY: str
    ALGORITHM: str = "HS256"
    ACCESS_TOKEN_EXPIRE_MINUTES: int = 30
    REFRESH_TOKEN_EXPIRE_DAYS: int = 7
    
    # Database
    DATABASE_URL: PostgresDsn
    DATABASE_POOL_SIZE: int = 20
    DATABASE_MAX_OVERFLOW: int = 0
    
    # Redis
    REDIS_URL: str
    REDIS_CACHE_TTL: int = 300
    
    # MinIO / S3
    MINIO_ENDPOINT: str
    MINIO_ACCESS_KEY: str
    MINIO_SECRET_KEY: str
    MINIO_BUCKET: str = "openbiobench"
    MINIO_SECURE: bool = False
    
    # CORS
    BACKEND_CORS_ORIGINS: List[str] = ["http://localhost:5173"]
    
    # File Upload
    MAX_UPLOAD_SIZE: int = 104857600  # 100MB
    
    # Scientific Tools
    AUTODOCK_VINA_IMAGE: str = "openbiobench/autodock-vina:1.2.3"
    RDKIT_IMAGE: str = "openbiobench/rdkit:2023.09"
    MGLTOOLS_IMAGE: str = "openbiobench/mgltools:1.5.7"
    GROMACS_IMAGE: str = "openbiobench/gromacs:2023.1"
    MAFFT_IMAGE: str = "biocontainers/mafft:v7.520"
    IQTREE_IMAGE: str = "biocontainers/iqtree:v2.2.0"
    
    # Job Management
    JOB_TIMEOUT_SECONDS: int = 3600
    MAX_CONCURRENT_JOBS: int = 10
    
    # External APIs
    UNIPROT_API_URL: str = "https://rest.uniprot.org"
    PDB_API_URL: str = "https://data.rcsb.org/rest/v1"
    ALPHAFOLD_API_URL: str = "https://alphafold.ebi.ac.uk/api"
    PUBCHEM_API_URL: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    CHEMBL_API_URL: str = "https://www.ebi.ac.uk/chembl/api/data"
    
    # Storage Paths (local development)
    STORAGE_ROOT: str = "./storage"
    UPLOAD_DIR: str = "./storage/uploads"
    RESULTS_DIR: str = "./storage/results"
    TEMP_DIR: str = "./storage/temp"


# Global settings instance
settings = Settings()
