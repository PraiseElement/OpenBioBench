"""
Job Queue Management Service.
Handles asynchronous job execution for computational tasks.
"""
import logging
import asyncio
from typing import Optional, Dict, Any
from datetime import datetime
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy import select
import uuid

from app.models.job import Job, JobStatus, JobType
from app.models.result import Result
from app.core.database import async_session_maker

logger = logging.getLogger(__name__)


class JobQueueService:
    """Service for managing computational job queue."""
    
    def __init__(self):
        self.running_jobs: Dict[str, asyncio.Task] = {}
    
    async def create_job(
        self,
        db: AsyncSession,
        user_id: uuid.UUID,
        project_id: uuid.UUID,
        job_type: JobType,
        parameters: Dict[str, Any],
        inputs: list[str]
    ) -> Job:
        """
        Create a new job and add to queue.
        
        Args:
            db: Database session
            user_id: User ID
            project_id: Project ID
            job_type: Type of job
            parameters: Job parameters
            inputs: Input entity IDs
            
        Returns:
            Created job
        """
        job = Job(
            user_id=user_id,
            project_id=project_id,
            job_type=job_type,
            status=JobStatus.PENDING,
            parameters=parameters,
            inputs=inputs,
            provenance={
                'created_by': str(user_id),
                'tool_version': '1.0.0',
                'parameters': parameters
            }
        )
        
        db.add(job)
        await db.commit()
        await db.refresh(job)
        
        logger.info(f"Created job {job.id} of type {job_type}")
        
        # Start job execution asynchronously
        task = asyncio.create_task(self._execute_job(str(job.id)))
        self.running_jobs[str(job.id)] = task
        
        return job
    
    async def _execute_job(self, job_id: str):
        """
        Execute a job asynchronously.
        
        Args:
            job_id: Job ID to execute
        """
        try:
            async with async_session_maker() as db:
                # Get job
                result = await db.execute(select(Job).where(Job.id == job_id))
                job = result.scalar_one_or_none()
                
                if not job:
                    logger.error(f"Job {job_id} not found")
                    return
                
                # Update status to running
                job.status = JobStatus.RUNNING
                job.started_at = datetime.utcnow()
                job.progress = 10
                await db.commit()
                
                logger.info(f"Starting execution of job {job_id} ({job.job_type})")
                
                # Execute based on job type
                try:
                    result_data = await self._run_job_logic(job, db)
                    
                    # Create result
                    result = Result(
                        job_id=job.id,
                        result_type=str(job.job_type),
                        summary=result_data.get('summary', {}),
                        data=result_data.get('data', {})
                    )
                    
                    db.add(result)
                    
                    # Update job as completed
                    job.status = JobStatus.COMPLETED
                    job.completed_at = datetime.utcnow()
                    job.progress = 100
                    job.compute_time_seconds = int(
                        (job.completed_at - job.started_at).total_seconds()
                    )
                    job.outputs = [str(result.id)]
                    
                    await db.commit()
                    
                    logger.info(f"Job {job_id} completed successfully")
                    
                except Exception as e:
                    # Job failed
                    job.status = JobStatus.FAILED
                    job.completed_at = datetime.utcnow()
                    job.error_message = str(e)
                    await db.commit()
                    
                    logger.error(f"Job {job_id} failed: {e}", exc_info=True)
                
        except Exception as e:
            logger.error(f"Error executing job {job_id}: {e}", exc_info=True)
        
        finally:
            # Remove from running jobs
            if job_id in self.running_jobs:
                del self.running_jobs[job_id]
    
    async def _run_job_logic(self, job: Job, db: AsyncSession) -> Dict:
        """
        Run job-specific logic.
        
        Args:
            job: Job to run
            db: Database session
            
        Returns:
            Result data dictionary
        """
        # Simulate processing time
        await asyncio.sleep(2)
        
        # Update progress
        job.progress = 50
        await db.commit()
        
        # Simple placeholder for different job types
        if job.job_type == JobType.LIGAND_PREPARATION:
            return {
                'summary': {'status': 'prepared', 'num_conformers': 10},
                'data': {'message': 'Ligand preparation completed'}
            }
        
        elif job.job_type == JobType.DOCKING:
            return {
                'summary': {
                    'status': 'completed',
                    'num_poses': 9,
                    'best_score': -8.5
                },
                'data': {
                    'poses': [
                        {'pose_id': 1, 'score': -8.5},
                        {'pose_id': 2, 'score': -8.2},
                    ]
                }
            }
        
        elif job.job_type == JobType.ADMET_PREDICTION:
            return {
                'summary': {'status': 'predicted'},
                'data': {'properties': {}}
            }
        
        else:
            return {
                'summary': {'status': 'completed'},
                'data': {}
            }
    
    async def get_job_status(self, db: AsyncSession, job_id: str) -> Optional[Job]:
        """
        Get job status.
        
        Args:
            db: Database session
            job_id: Job ID
            
        Returns:
            Job object or None
        """
        result = await db.execute(select(Job).where(Job.id == job_id))
        return result.scalar_one_or_none()
    
    async def cancel_job(self, db: AsyncSession, job_id: str) -> bool:
        """
        Cancel a running job.
        
        Args:
            db: Database session
            job_id: Job ID
            
        Returns:
            Success status
        """
        # Get job
        result = await db.execute(select(Job).where(Job.id == job_id))
        job = result.scalar_one_or_none()
        
        if not job:
            return False
        
        # Cancel if running
        if job.status in [JobStatus.PENDING, JobStatus.RUNNING]:
            if job_id in self.running_jobs:
                self.running_jobs[job_id].cancel()
                del self.running_jobs[job_id]
            
            job.status = JobStatus.CANCELLED
            job.completed_at = datetime.utcnow()
            await db.commit()
            
            logger.info(f"Job {job_id} cancelled")
            return True
        
        return False


# Global job queue instance
job_queue = JobQueueService()
