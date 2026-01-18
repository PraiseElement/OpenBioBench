/**
 * TypeScript type definitions for OpenBioBench.
 * Mirrors backend Pydantic models.
 */

// User types
export interface User {
  id: string;
  email: string;
  institution?: string;
  role: 'admin' | 'researcher' | 'guest';
  compute_quota_hours: number;
  storage_quota_gb: number;
  is_active: boolean;
  created_at: string;
  last_login?: string;
}

export interface TokenResponse {
  access_token: string;
  refresh_token: string;
  token_type: string;
  expires_in: number;
}

// Ligand types
export interface Ligand {
  id: string;
  project_id: string;
  source: string;
  smiles: string;
  inchi?: string;
  inchi_key?: string;
  name?: string;
  mol_file_url?: string;
  properties: MolecularProperties;
  metadata: Record<string, any>;
  created_at: string;
}

export interface MolecularProperties {
  molecular_weight: number;
  logp: number;
  tpsa: number;
  hbd: number;
  hba: number;
  rotatable_bonds: number;
  aromatic_rings: number;
  lipinski_violations: number;
  molecular_formula?: string;
}

// Builder types
export interface ConformerInfo {
  conformer_id: number;
  energy: number;
  rmsd?: number;
  file_url: string;
}

export interface Generate3DResponse {
  ligand_id: string;
  smiles: string;
  canonical_smiles: string;
  conformers: ConformerInfo[];
  lowest_energy_conformer_id: number;
  properties: MolecularProperties;
}

// ADMET types
export interface PropertyPrediction {
  value: number;
  unit?: string;
  confidence: 'low' | 'medium' | 'high';
  interpretation?: string;
}

export interface PhysicochemicalProperties {
  molecular_weight: number;
  logp: number;
  logd?: number;
  tpsa: number;
  hbd: number;
  hba: number;
  rotatable_bonds: number;
  aromatic_rings: number;
  lipinski_violations: number;
}

export interface PharmacokineticProperties {
  caco2_permeability: PropertyPrediction;
  bbb_penetration: PropertyPrediction;
  pgp_substrate?: PropertyPrediction;
  hia?: PropertyPrediction;
}

export interface ToxicityProperties {
  herg_ic50: PropertyPrediction;
  ames_mutagenicity: {
    prediction: string;
    confidence: number;
  };
  hepatotoxicity?: any;
  ld50?: PropertyPrediction;
}

export interface ADMETResponse {
  ligand_id?: string;
  smiles: string;
  canonical_smiles: string;
  properties: PhysicochemicalProperties;
  pharmacokinetics: PharmacokineticProperties;
  toxicity: ToxicityProperties;
  warnings: string[];
}

// Job types
export interface Job {
  id: string;
  user_id: string;
  project_id: string;
  job_type: string;
  status: 'pending' | 'queued' | 'running' | 'completed' | 'failed' | 'cancelled';
  parameters: Record<string, any>;
  inputs: string[];
  outputs: string[];
  created_at: string;
  started_at?: string;
  completed_at?: string;
  compute_time_seconds?: number;
  error_message?: string;
  provenance: Record<string, any>;
  progress: number;
}

// Project types
export interface Project {
  id: string;
  owner_id: string;
  name: string;
  description?: string;
  tags: string[];
  created_at: string;
  updated_at: string;
}
