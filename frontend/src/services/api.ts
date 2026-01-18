/**
 * API client for OpenBioBench backend.
 * Provides typed methods for all API endpoints.
 */
import axios, { type AxiosInstance, type AxiosError } from 'axios';
import type {
  User,
  TokenResponse,
  Ligand,
  Generate3DResponse,
  ADMETResponse,
  Job,
  Project,
} from '../types';

const API_BASE_URL = import.meta.env.VITE_API_URL 
  ? `${import.meta.env.VITE_API_URL}/api` 
  : 'http://localhost:8000/api';

class ApiClient {
  private client: AxiosInstance;
  private accessToken: string | null = null;

  constructor() {
    this.client = axios.create({
      baseURL: API_BASE_URL,
      headers: {
        'Content-Type': 'application/json',
      },
    });

    // Request interceptor to add auth token
    this.client.interceptors.request.use((config) => {
      if (this.accessToken) {
        config.headers.Authorization = `Bearer ${this.accessToken}`;
      }
      return config;
    });

    // Response interceptor for error handling
    this.client.interceptors.response.use(
      (response) => response,
      async (error: AxiosError) => {
        if (error.response?.status === 401) {
          // Token expired, try to refresh
          // For now, just clear token and redirect to login
          this.setToken(null);
          window.location.href = '/login';
        }
        return Promise.reject(error);
      }
    );

    // Load token from localStorage
    const storedToken = localStorage.getItem('access_token');
    if (storedToken) {
      this.accessToken = storedToken;
    }
  }

  setToken(token: string | null) {
    this.accessToken = token;
    if (token) {
      localStorage.setItem('access_token', token);
    } else {
      localStorage.removeItem('access_token');
    }
  }

  // Authentication
  async register(email: string, password: string, institution?: string): Promise<User> {
    const { data } = await this.client.post<User>('/v1/auth/register', {
      email,
      password,
      institution,
    });
    return data;
  }

  async login(email: string, password: string): Promise<TokenResponse> {
    const formData = new FormData();
    formData.append('username', email);
    formData.append('password', password);

    const { data } = await this.client.post<TokenResponse>('/v1/auth/login', formData, {
      headers: {
        'Content-Type': 'application/x-www-form-urlencoded',
      },
    });

    this.setToken(data.access_token);
    localStorage.setItem('refresh_token', data.refresh_token);

    return data;
  }

  async getCurrentUser(): Promise<User> {
    const { data } = await this.client.get<User>('/v1/auth/me');
    return data;
  }

  async logout() {
    this.setToken(null);
    localStorage.removeItem('refresh_token');
  }

  // Builder
  async validateSMILES(smiles: string) {
    const { data } = await this.client.post('/v1/builder/validate', { smiles });
    return data;
  }

  async generate3D(params: {
    smiles: string;
    project_id: string;
    method?: string;
    num_conformers?: number;
    minimize?: boolean;
    force_field?: string;
    name?: string;
  }): Promise<Generate3DResponse> {
    const { data } = await this.client.post<Generate3DResponse>('/v1/builder/2d/to_3d', params);
    return data;
  }

  // ADMET
  async predictADMET(params: { smiles?: string; ligand_id?: string }): Promise<ADMETResponse> {
    const { data } = await this.client.post<ADMETResponse>('/v1/admet/predict', params);
    return data;
  }

  // Projects (placeholder)
  async getProjects(): Promise<Project[]> {
    // Implementation needed
    return [];
  }

  async createProject(name: string, description?: string, tags?: string[]): Promise<Project> {
    // Implementation needed
    return {} as Project;
  }

  // Jobs (placeholder)
  async getJobs(): Promise<Job[]> {
    // Implementation needed
    return [];
  }

  async getJob(jobId: string): Promise<Job> {
    // Implementation needed
    return {} as Job;
  }
}

export const apiClient = new ApiClient();
