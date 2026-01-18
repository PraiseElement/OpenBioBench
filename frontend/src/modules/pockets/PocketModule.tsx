/**
 * Pocket Detection Module.
 * Binding site identification with database fetch.
 */
import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import { Target, ArrowLeft, Loader, Download, Upload, Search } from 'lucide-react';

export const PocketModule: React.FC = () => {
  const [pdbContent, setPdbContent] = useState('');
  const [pdbId, setPdbId] = useState('');
  const [pockets, setPockets] = useState<any[]>([]);
  const [isLoading, setIsLoading] = useState(false);
  const [isFetching, setIsFetching] = useState(false);
  const [error, setError] = useState('');

  const fetchFromPDB = async () => {
    if (!pdbId.trim()) return;
    setIsFetching(true);
    setError('');
    
    try {
      const response = await fetch(`http://localhost:8000/api/v1/database/pdb/${pdbId}`, {
        headers: { 'Authorization': `Bearer ${localStorage.getItem('access_token')}` },
      });
      
      if (!response.ok) throw new Error(`PDB ID ${pdbId} not found`);
      
      const data = await response.json();
      setPdbContent(data.pdb_content);
    } catch (err: any) {
      setError(err.message || 'Failed to fetch from PDB');
    } finally {
      setIsFetching(false);
    }
  };

  const fetchFromAlphaFold = async (uniprotId: string) => {
    setIsFetching(true);
    setError('');
    
    try {
      const response = await fetch(`http://localhost:8000/api/v1/database/alphafold/${uniprotId}`, {
        headers: { 'Authorization': `Bearer ${localStorage.getItem('access_token')}` },
      });
      
      if (!response.ok) throw new Error('AlphaFold structure not found');
      
      const data = await response.json();
      setPdbContent(data.pdb_content);
    } catch (err: any) {
      setError(err.message || 'Failed to fetch from AlphaFold');
    } finally {
      setIsFetching(false);
    }
  };

  const handleFileUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;
    
    const content = await file.text();
    setPdbContent(content);
  };

  const handleDetect = async () => {
    if (!pdbContent.trim()) return;

    setIsLoading(true);
    setError('');
    setPockets([]);

    try {
      const response = await fetch('http://localhost:8000/api/v1/pockets/detect', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${localStorage.getItem('access_token')}`,
        },
        body: JSON.stringify({ pdb_content: pdbContent }),
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Detection failed');
      }

      const data = await response.json();
      setPockets(data.pockets || []);
    } catch (err: any) {
      setError(err.message || 'Detection failed');
    } finally {
      setIsLoading(false);
    }
  };

  const getDruggabilityColor = (score: number) => {
    if (score >= 0.7) return 'bg-green-500';
    if (score >= 0.5) return 'bg-yellow-500';
    if (score >= 0.3) return 'bg-orange-500';
    return 'bg-red-500';
  };

  const getDruggabilityLabel = (score: number) => {
    if (score >= 0.7) return 'High';
    if (score >= 0.5) return 'Moderate';
    if (score >= 0.3) return 'Low';
    return 'Very Low';
  };

  return (
    <div className="min-h-screen bg-slate-50 dark:bg-slate-900">
      {/* Header */}
      <header className="bg-white dark:bg-slate-800 border-b dark:border-slate-700">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-4">
          <div className="flex items-center space-x-4">
            <Link to="/" className="text-slate-600 hover:text-slate-900">
              <ArrowLeft className="h-6 w-6" />
            </Link>
            <div className="flex items-center space-x-3">
              <Target className="h-8 w-8 text-orange-600" />
              <div>
                <h1 className="text-2xl font-bold text-slate-900 dark:text-slate-100">Pocket Detection</h1>
                <p className="text-sm text-slate-600 dark:text-slate-400">Binding Site Identification & Druggability</p>
              </div>
            </div>
          </div>
        </div>
      </header>

      <main className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
          {/* Input Panel */}
          <div className="card">
            <h2 className="text-xl font-semibold mb-4">Protein Structure</h2>

            {/* Fetch from PDB */}
            <div className="mb-4">
              <label className="block text-sm font-medium text-slate-700 mb-1">
                Fetch from RCSB PDB
              </label>
              <div className="flex space-x-2">
                <input
                  type="text"
                  className="input flex-1"
                  value={pdbId}
                  onChange={(e) => setPdbId(e.target.value.toUpperCase())}
                  placeholder="e.g., 1HSG, 4HHB, 6LU7"
                />
                <button
                  onClick={fetchFromPDB}
                  className="btn btn-secondary flex items-center space-x-2"
                  disabled={isFetching || !pdbId.trim()}
                >
                  <Search className="h-4 w-4" />
                  <span>Fetch</span>
                </button>
              </div>
            </div>

            {/* Well-known examples */}
            <div className="mb-4">
              <p className="text-sm text-slate-500 mb-2">Quick examples:</p>
              <div className="flex gap-2 flex-wrap">
                <button
                  onClick={() => { setPdbId('1HSG'); }}
                  className="text-xs px-3 py-1 bg-slate-100 rounded-full hover:bg-slate-200"
                >
                  HIV Protease (1HSG)
                </button>
                <button
                  onClick={() => { setPdbId('4HHB'); }}
                  className="text-xs px-3 py-1 bg-slate-100 rounded-full hover:bg-slate-200"
                >
                  Hemoglobin (4HHB)
                </button>
                <button
                  onClick={() => { setPdbId('6LU7'); }}
                  className="text-xs px-3 py-1 bg-slate-100 rounded-full hover:bg-slate-200"
                >
                  SARS-CoV-2 (6LU7)
                </button>
              </div>
            </div>

            {/* File upload */}
            <div className="mb-4">
              <label className="block text-sm font-medium text-slate-700 mb-1">
                Or Upload PDB File
              </label>
              <input
                type="file"
                onChange={handleFileUpload}
                accept=".pdb"
                className="block w-full text-sm text-slate-500 file:mr-4 file:py-2 file:px-4 file:rounded-lg file:border-0 file:text-sm file:font-medium file:bg-orange-50 file:text-orange-700 hover:file:bg-orange-100"
              />
            </div>

            {pdbContent && (
              <div className="mb-4 p-3 bg-green-50 border border-green-200 rounded-lg text-sm text-green-700">
                ✓ Structure loaded ({pdbContent.split('\n').filter(l => l.startsWith('ATOM')).length} atoms)
              </div>
            )}

            <button
              onClick={handleDetect}
              className="w-full btn btn-primary"
              disabled={isLoading || !pdbContent.trim()}
            >
              {isLoading ? (
                <span className="flex items-center justify-center space-x-2">
                  <Loader className="h-4 w-4 animate-spin" />
                  <span>Detecting Pockets...</span>
                </span>
              ) : (
                'Detect Binding Pockets'
              )}
            </button>

            {error && (
              <div className="mt-4 bg-red-50 border border-red-200 text-red-700 px-4 py-3 rounded-lg text-sm">
                {error}
              </div>
            )}
          </div>

          {/* Results Panel */}
          <div className="card">
            <h2 className="text-xl font-semibold mb-4">Detected Pockets</h2>

            {pockets.length === 0 && !error && (
              <div className="text-center py-12 text-slate-500">
                <Target className="h-16 w-16 mx-auto mb-4 text-slate-300" />
                <p>Fetch or upload a protein structure to detect binding pockets</p>
              </div>
            )}

            {pockets.length > 0 && (
              <div className="space-y-4">
                <div className="bg-orange-50 border border-orange-200 rounded-lg p-4">
                  <p className="font-medium text-orange-900">
                    {pockets.length} pocket(s) detected
                  </p>
                  <p className="text-sm text-orange-700">
                    Ranked by druggability score
                  </p>
                </div>

                <div className="space-y-3">
                  {pockets.map((pocket, idx) => (
                    <div key={idx} className="p-4 border rounded-lg hover:shadow-md transition-shadow">
                      <div className="flex items-center justify-between mb-2">
                        <h3 className="font-semibold">Pocket {pocket.pocket_id}</h3>
                        <div className="flex items-center space-x-2">
                          <div className={`w-3 h-3 rounded-full ${getDruggabilityColor(pocket.druggability_score)}`}></div>
                          <span className="text-sm font-medium">
                            {getDruggabilityLabel(pocket.druggability_score)}
                          </span>
                        </div>
                      </div>
                      
                      <div className="grid grid-cols-2 gap-2 text-sm text-slate-600">
                        <div>Druggability: {(pocket.druggability_score * 100).toFixed(0)}%</div>
                        <div>Volume: {pocket.volume?.toFixed(0)} Å³</div>
                        <div>Surface: {pocket.surface_area?.toFixed(0)} Å²</div>
                        <div>Depth: {pocket.depth?.toFixed(1)} Å</div>
                        <div>Hydrophobicity: {(pocket.hydrophobicity_ratio * 100).toFixed(0)}%</div>
                        <div>Center: [{pocket.center?.map((c: number) => c.toFixed(1)).join(', ')}]</div>
                      </div>

                      {pocket.residues && pocket.residues.length > 0 && (
                        <div className="mt-2 text-xs text-slate-500">
                          Residues: {pocket.residues.join(', ')}
                        </div>
                      )}
                    </div>
                  ))}
                </div>

                <button className="btn btn-secondary flex items-center space-x-2">
                  <Download className="h-4 w-4" />
                  <span>Export Results</span>
                </button>
              </div>
            )}
          </div>
        </div>
      </main>
    </div>
  );
};
