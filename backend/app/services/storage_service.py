"""
Storage service for handling file uploads to MinIO/S3.
Manages molecular structure files, results, and user data.
"""
import os
import logging
import uuid
from typing import Optional, BinaryIO
from pathlib import Path
from minio import Minio
from minio.error import S3Error
from app.core.config import settings

logger = logging.getLogger(__name__)


class StorageService:
    """Service for file storage operations."""
    
    def __init__(self):
        """Initialize MinIO client."""
        self.client = None
        self.bucket = settings.MINIO_BUCKET
        
        # Initialize MinIO client if configured
        if settings.MINIO_ENDPOINT:
            try:
                self.client = Minio(
                    settings.MINIO_ENDPOINT,
                    access_key=settings.MINIO_ACCESS_KEY,
                    secret_key=settings.MINIO_SECRET_KEY,
                    secure=settings.MINIO_SECURE
                )
                
                # Ensure bucket exists
                if not self.client.bucket_exists(self.bucket):
                    self.client.make_bucket(self.bucket)
                    logger.info(f"Created MinIO bucket: {self.bucket}")
                    
            except Exception as e:
                logger.warning(f"MinIO not available, using local storage: {e}")
                self.client = None
        
        # Ensure local storage directories exist
        for dir_path in [settings.UPLOAD_DIR, settings.RESULTS_DIR, settings.TEMP_DIR]:
            Path(dir_path).mkdir(parents=True, exist_ok=True)
    
    def upload_file(
        self,
        file_data: bytes,
        filename: str,
        content_type: str = "application/octet-stream",
        folder: str = "uploads"
    ) -> str:
        """
        Upload file to storage.
        
        Args:
            file_data: File binary data
            filename: Original filename
            content_type: MIME type
            folder: Folder/prefix in storage
            
        Returns:
            Storage URL/path for the uploaded file
        """
        # Generate unique filename
        file_ext = Path(filename).suffix
        unique_name = f"{uuid.uuid4()}{file_ext}"
        object_name = f"{folder}/{unique_name}"
        
        if self.client:
            # Upload to MinIO
            try:
                from io import BytesIO
                self.client.put_object(
                    self.bucket,
                    object_name,
                    BytesIO(file_data),
                    length=len(file_data),
                    content_type=content_type
                )
                return f"s3://{self.bucket}/{object_name}"
                
            except S3Error as e:
                logger.error(f"MinIO upload failed: {e}")
                # Fall through to local storage
        
        # Use local storage
        local_path = Path(settings.STORAGE_ROOT) / object_name
        local_path.parent.mkdir(parents=True, exist_ok=True)
        local_path.write_bytes(file_data)
        
        return str(local_path)
    
    def upload_text(
        self,
        text: str,
        filename: str,
        folder: str = "results"
    ) -> str:
        """Upload text file to storage."""
        return self.upload_file(
            text.encode('utf-8'),
            filename,
            "text/plain",
            folder
        )
    
    def download_file(self, file_url: str) -> Optional[bytes]:
        """
        Download file from storage.
        
        Args:
            file_url: Storage URL or local path
            
        Returns:
            File binary data or None if not found
        """
        if file_url.startswith("s3://") and self.client:
            # Download from MinIO
            try:
                bucket_name, object_name = file_url.replace("s3://", "").split("/", 1)
                response = self.client.get_object(bucket_name, object_name)
                data = response.read()
                response.close()
                response.release_conn()
                return data
                
            except S3Error as e:
                logger.error(f"MinIO download failed: {e}")
                return None
        else:
            # Read from local storage
            try:
                return Path(file_url).read_bytes()
            except Exception as e:
                logger.error(f"Local file read failed: {e}")
                return None
    
    def delete_file(self, file_url: str) -> bool:
        """Delete file from storage."""
        if file_url.startswith("s3://") and self.client:
            try:
                bucket_name, object_name = file_url.replace("s3://", "").split("/", 1)
                self.client.remove_object(bucket_name, object_name)
                return True
            except S3Error as e:
                logger.error(f"MinIO delete failed: {e}")
                return False
        else:
            try:
                Path(file_url).unlink(missing_ok=True)
                return True
            except Exception as e:
                logger.error(f"Local file delete failed: {e}")
                return False
    
    def get_presigned_url(self, file_url: str, expires: int = 3600) -> Optional[str]:
        """
        Generate presigned URL for file download.
        
        Args:
            file_url: Storage URL
            expires: URL expiration in seconds
            
        Returns:
            Presigned URL or original path
        """
        if file_url.startswith("s3://") and self.client:
            try:
                bucket_name, object_name = file_url.replace("s3://", "").split("/", 1)
                url = self.client.presigned_get_object(
                    bucket_name,
                    object_name,
                    expires=expires
                )
                return url
            except S3Error as e:
                logger.error(f"Presigned URL generation failed: {e}")
                return file_url
        else:
            # For local files, return the path directly
            return file_url


# Global storage instance
storage = StorageService()
