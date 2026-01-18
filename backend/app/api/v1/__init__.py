"""API v1 router initialization."""
from fastapi import APIRouter
from app.api.v1 import (
    auth, builder, admet, projects, sequences, structures, jobs,
    alignments, pockets, docking, phylogeny, database, files, preparation
)

api_router = APIRouter(prefix="/v1")

# Include all endpoint routers
api_router.include_router(auth.router)
api_router.include_router(builder.router)
api_router.include_router(admet.router)
api_router.include_router(projects.router)
api_router.include_router(sequences.router)
api_router.include_router(structures.router)
api_router.include_router(jobs.router)
api_router.include_router(alignments.router)
api_router.include_router(pockets.router)
api_router.include_router(docking.router)
api_router.include_router(phylogeny.router)
api_router.include_router(database.router)
api_router.include_router(files.router)
api_router.include_router(preparation.router)

__all__ = ["api_router"]


