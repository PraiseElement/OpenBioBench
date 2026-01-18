/**
 * API configuration and helper utilities.
 * Centralized API URL management.
 */

// Get the base API URL from environment variable or fallback to localhost
export const API_BASE_URL = import.meta.env.VITE_API_URL 
  ? `${import.meta.env.VITE_API_URL}/api` 
  : 'http://localhost:8000/api';

// Helper function to get full API URL
export const getApiUrl = (path: string): string => {
  // Remove leading slash if present
  const cleanPath = path.startsWith('/') ? path.slice(1) : path;
  return `${API_BASE_URL}/${cleanPath}`;
};

// Export for backward compatibility
export default API_BASE_URL;
