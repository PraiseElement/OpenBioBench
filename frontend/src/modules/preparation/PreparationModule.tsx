/**
 * Structure Preparation Module.
 * Clean and prepare proteins and ligands for docking.
 */
import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import { Wrench, ArrowLeft, Loader, Download, Upload, Search } from 'lucide-react';

export const PreparationModule: React.FC = () => {
  const [activeTab, setActiveTab] = useState<'protein' | 'ligand'>('protein');
  
  // Protein state
  const [proteinPdb, setProteinPdb] = useState('');
  const [pdbId, setPdbId] = useState('');
  const [removeWater, setRemoveWater] = useState(true);
  const [removeHeteroatoms, setRemoveHeteroatoms] = useState(false);
  const [cleanedProtein, setCleanedProtein] = useState<any>(null);
  
  // Ligand state
  const [ligandSmiles, setLigandSmiles] = useState('');
  const [ligandFile, setLigandFile] = useState<File | null>(null);
  const [preparedLigand, setPreparedLigand] = useState<any>(null);
  
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState('');

  const fetchFromPDB = async () => {
    if (!pdbId.trim()) return;
    setIsLoading(true);
    setError('');
    
    try {
      const response = await fetch(`http://localhost:8000/api/v1/database/pdb/${pdbId}`, {
        headers: { 'Authorization': `Bearer ${localStorage.getItem('access_token')}` },
      });
      
      if (!response.ok) throw new Error('PDB ID not found');
      
      const data = await response.json();
      setProteinPdb(data.pdb_content);
    } catch (err: any) {
      setError(err.message || 'Failed to fetch from PDB');
    } finally {
      setIsLoading(false);
    }
  };

  const cleanProtein = async () => {
    if (!proteinPdb.trim()) return;
    setIsLoading(true);
    setError('');
    setCleanedProtein(null);
    
    try {
      const response = await fetch('http://localhost:8000/api/v1/preparation/protein/clean', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${localStorage.getItem('access_token')}`,
        },
        body: JSON.stringify({
          pdb_content: proteinPdb,
          remove_water: removeWater,
          remove_heteroatoms: removeHeteroatoms,
        }),
      });
      
      if (!response.ok) throw new Error('Protein cleaning failed');
      
      const data = await response.json();
      setCleanedProtein(data);
    } catch (err: any) {
      setError(err.message || 'Cleaning failed');
    } finally {
      setIsLoading(false);
    }
  };

  const prepareLigand = async () => {
    if (!ligandSmiles.trim()) return;
    setIsLoading(true);
    setError('');
    setPreparedLigand(null);
    
    try {
      const response = await fetch('http://localhost:8000/api/v1/preparation/ligand/from_smiles', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${localStorage.getItem('access_token')}`,
        },
        body: JSON.stringify({
          smiles: ligandSmiles,
          optimize: true,
        }),
      });
      
      if (!response.ok) throw new Error('Ligand preparation failed');
      
      const data = await response.json();
      setPreparedLigand(data);
    } catch (err: any) {
      setError(err.message || 'Preparation failed');
    } finally {
      setIsLoading(false);
    }
  };

  const downloadFile = (content: string, filename: string) => {
    const blob = new Blob([content], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
  };

  return (
    <div className="min-h-screen bg-slate-50">
      {/* Header */}
      <header className="bg-white border-b">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-4">
          <div className="flex items-center space-x-4">
            <Link to="/" className="text-slate-600 hover:text-slate-900">
              <ArrowLeft className="h-6 w-6" />
            </Link>
            <div className="flex items-center space-x-3">
              <Wrench className="h-8 w-8 text-amber-600" />
              <div>
                <h1 className="text-2xl font-bold text-slate-900">Structure Preparation</h1>
                <p className="text-sm text-slate-600">Clean & Prepare Proteins and Ligands</p>
              </div>
            </div>
          </div>
        </div>
      </header>

      <main className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        {/* Tabs */}
        <div className="flex space-x-4 mb-6">
          <button
            onClick={() => setActiveTab('protein')}
            className={`px-6 py-2 rounded-lg font-medium transition-colors ${
              activeTab === 'protein'
                ? 'bg-amber-500 text-white'
                : 'bg-white text-slate-600 hover:bg-slate-100'
            }`}
          >
            Protein Preparation
          </button>
          <button
            onClick={() => setActiveTab('ligand')}
            className={`px-6 py-2 rounded-lg font-medium transition-colors ${
              activeTab === 'ligand'
                ? 'bg-amber-500 text-white'
                : 'bg-white text-slate-600 hover:bg-slate-100'
            }`}
          >
            Ligand Preparation
          </button>
        </div>

        {error && (
          <div className="mb-6 bg-red-50 border border-red-200 text-red-700 px-4 py-3 rounded-lg">
            {error}
          </div>
        )}

        {activeTab === 'protein' && (
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
            {/* Input */}
            <div className="card">
              <h2 className="text-xl font-semibold mb-4">Protein Structure</h2>
              
              {/* Fetch from PDB */}
              <div className="mb-4">
                <label className="block text-sm font-medium text-slate-700 mb-1">
                  Fetch from PDB
                </label>
                <div className="flex space-x-2">
                  <input
                    type="text"
                    className="input flex-1"
                    value={pdbId}
                    onChange={(e) => setPdbId(e.target.value.toUpperCase())}
                    placeholder="e.g., 1HSG, 4HHB"
                  />
                  <button
                    onClick={fetchFromPDB}
                    className="btn btn-secondary flex items-center space-x-2"
                    disabled={isLoading || !pdbId.trim()}
                  >
                    <Search className="h-4 w-4" />
                    <span>Fetch</span>
                  </button>
                </div>
              </div>

              <div className="text-center text-slate-500 text-sm my-4">— or paste PDB content —</div>

              <textarea
                className="input font-mono text-xs"
                rows={10}
                value={proteinPdb}
                onChange={(e) => setProteinPdb(e.target.value)}
                placeholder="Paste PDB content here..."
              />

              {/* Options */}
              <div className="mt-4 space-y-2">
                <label className="flex items-center space-x-2">
                  <input
                    type="checkbox"
                    checked={removeWater}
                    onChange={(e) => setRemoveWater(e.target.checked)}
                    className="rounded"
                  />
                  <span className="text-sm">Remove water molecules (HOH)</span>
                </label>
                <label className="flex items-center space-x-2">
                  <input
                    type="checkbox"
                    checked={removeHeteroatoms}
                    onChange={(e) => setRemoveHeteroatoms(e.target.checked)}
                    className="rounded"
                  />
                  <span className="text-sm">Remove all heteroatoms (ligands, ions)</span>
                </label>
              </div>

              <button
                onClick={cleanProtein}
                className="w-full btn btn-primary mt-4"
                disabled={isLoading || !proteinPdb.trim()}
              >
                {isLoading ? (
                  <span className="flex items-center justify-center space-x-2">
                    <Loader className="h-4 w-4 animate-spin" />
                    <span>Processing...</span>
                  </span>
                ) : (
                  'Clean & Prepare Protein'
                )}
              </button>
            </div>

            {/* Output */}
            <div className="card">
              <h2 className="text-xl font-semibold mb-4">Prepared Structure</h2>
              
              {!cleanedProtein && (
                <div className="text-center py-12 text-slate-500">
                  <Wrench className="h-16 w-16 mx-auto mb-4 text-slate-300" />
                  <p>Fetch or paste a protein structure to prepare it</p>
                </div>
              )}

              {cleanedProtein && (
                <div className="space-y-4">
                  <div className="bg-green-50 border border-green-200 rounded-lg p-4">
                    <h3 className="font-semibold text-green-900 mb-2">Cleaning Summary</h3>
                    <div className="grid grid-cols-2 gap-2 text-sm">
                      <div>Waters removed: {cleanedProtein.stats?.removed_waters || 0}</div>
                      <div>Heteroatoms removed: {cleanedProtein.stats?.removed_heteroatoms || 0}</div>
                      <div>Atoms kept: {cleanedProtein.stats?.kept_atoms || 0}</div>
                    </div>
                  </div>

                  <button
                    onClick={() => downloadFile(cleanedProtein.pdb_content, 'prepared_protein.pdb')}
                    className="btn btn-primary flex items-center space-x-2"
                  >
                    <Download className="h-4 w-4" />
                    <span>Download Prepared PDB</span>
                  </button>
                </div>
              )}
            </div>
          </div>
        )}

        {activeTab === 'ligand' && (
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
            {/* Input */}
            <div className="card">
              <h2 className="text-xl font-semibold mb-4">Ligand Input</h2>
              
              <div className="space-y-4">
                <div>
                  <label className="block text-sm font-medium text-slate-700 mb-1">
                    SMILES String
                  </label>
                  <input
                    type="text"
                    className="input font-mono"
                    value={ligandSmiles}
                    onChange={(e) => setLigandSmiles(e.target.value)}
                    placeholder="e.g., CC(=O)Oc1ccccc1C(=O)O (Aspirin)"
                  />
                </div>

                <div className="flex gap-2 flex-wrap">
                  <button
                    onClick={() => setLigandSmiles('CC(=O)Oc1ccccc1C(=O)O')}
                    className="text-xs px-3 py-1 bg-slate-100 rounded-full hover:bg-slate-200"
                  >
                    Aspirin
                  </button>
                  <button
                    onClick={() => setLigandSmiles('CN1C=NC2=C1C(=O)N(C(=O)N2C)C')}
                    className="text-xs px-3 py-1 bg-slate-100 rounded-full hover:bg-slate-200"
                  >
                    Caffeine
                  </button>
                  <button
                    onClick={() => setLigandSmiles('CC(C)Cc1ccc(C(C)C(=O)O)cc1')}
                    className="text-xs px-3 py-1 bg-slate-100 rounded-full hover:bg-slate-200"
                  >
                    Ibuprofen
                  </button>
                </div>

                <button
                  onClick={prepareLigand}
                  className="w-full btn btn-primary"
                  disabled={isLoading || !ligandSmiles.trim()}
                >
                  {isLoading ? (
                    <span className="flex items-center justify-center space-x-2">
                      <Loader className="h-4 w-4 animate-spin" />
                      <span>Preparing...</span>
                    </span>
                  ) : (
                    'Generate 3D & Optimize'
                  )}
                </button>
              </div>
            </div>

            {/* Output */}
            <div className="card">
              <h2 className="text-xl font-semibold mb-4">Prepared Ligand</h2>
              
              {!preparedLigand && (
                <div className="text-center py-12 text-slate-500">
                  <Wrench className="h-16 w-16 mx-auto mb-4 text-slate-300" />
                  <p>Enter a SMILES string to prepare a ligand</p>
                </div>
              )}

              {preparedLigand && (
                <div className="space-y-4">
                  <div className="bg-amber-50 border border-amber-200 rounded-lg p-4">
                    <h3 className="font-semibold text-amber-900 mb-2">Ligand Properties</h3>
                    <div className="grid grid-cols-2 gap-2 text-sm">
                      <div>MW: {preparedLigand.properties?.molecular_weight}</div>
                      <div>LogP: {preparedLigand.properties?.logP}</div>
                      <div>HBD: {preparedLigand.properties?.hbd}</div>
                      <div>HBA: {preparedLigand.properties?.hba}</div>
                      <div>Rotatable: {preparedLigand.properties?.rotatable_bonds}</div>
                      <div>Atoms: {preparedLigand.num_atoms}</div>
                    </div>
                  </div>

                  <div className="flex space-x-2">
                    <button
                      onClick={() => downloadFile(preparedLigand.pdb_content, 'prepared_ligand.pdb')}
                      className="btn btn-primary flex items-center space-x-2"
                    >
                      <Download className="h-4 w-4" />
                      <span>Download PDB</span>
                    </button>
                    <button
                      onClick={() => downloadFile(preparedLigand.mol_content, 'prepared_ligand.mol')}
                      className="btn btn-secondary flex items-center space-x-2"
                    >
                      <Download className="h-4 w-4" />
                      <span>Download MOL</span>
                    </button>
                  </div>
                </div>
              )}
            </div>
          </div>
        )}
      </main>
    </div>
  );
};
