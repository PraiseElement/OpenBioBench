/**
 * Molecular Docking Module.
 * Protein-ligand docking with database fetch and file upload.
 */
import React, { useState, useRef } from 'react';
import { Link } from 'react-router-dom';
import { Database, ArrowLeft, Loader, Upload, Search } from 'lucide-react';

export const DockingModule: React.FC = () => {
  const [proteinPdb, setProteinPdb] = useState('');
  const [ligandSmiles, setLigandSmiles] = useState('');
  const [boxCenter, setBoxCenter] = useState({ x: 0, y: 0, z: 0 });
  const [boxSize, setBoxSize] = useState({ x: 20, y: 20, z: 20 });
  const [poses, setPoses] = useState<any[]>([]);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState('');
  
  // Database/file state
  const [pdbId, setPdbId] = useState('');
  const [pubchemCid, setPubchemCid] = useState('');
  const [isFetching, setIsFetching] = useState(false);
  
  const proteinFileRef = useRef<HTMLInputElement>(null);
  const ligandFileRef = useRef<HTMLInputElement>(null);

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
      setProteinPdb(data.pdb_content);
      
      // Try to auto-calculate box if structure is loaded
      if (data.pdb_content) {
        calculateDockingBox(data.pdb_content);
      }
    } catch (err: any) {
      setError(err.message || 'Failed to fetch from PDB');
    } finally {
      setIsFetching(false);
    }
  };

  const fetchFromPubChem = async () => {
    if (!pubchemCid.trim()) return;
    setIsFetching(true);
    setError('');
    
    try {
      const response = await fetch(`http://localhost:8000/api/v1/database/pubchem/${pubchemCid}`, {
        headers: { 'Authorization': `Bearer ${localStorage.getItem('access_token')}` },
      });
      
      if (!response.ok) throw new Error(`PubChem CID ${pubchemCid} not found`);
      
      const data = await response.json();
      setLigandSmiles(data.smiles);
    } catch (err: any) {
      setError(err.message || 'Failed to fetch from PubChem');
    } finally {
      setIsFetching(false);
    }
  };

  const handleProteinFileUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;
    
    const content = await file.text();
    setProteinPdb(content);
    
    if (content.includes('ATOM')) {
      calculateDockingBox(content);
    }
  };

  const handleLigandFileUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    if (!file) return;
    
    const content = await file.text();
    
    // Try to prepare ligand and get SMILES
    try {
      const response = await fetch('http://localhost:8000/api/v1/preparation/ligand/from_sdf', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${localStorage.getItem('access_token')}`,
        },
        body: JSON.stringify({ sdf_content: content }),
      });
      
      if (response.ok) {
        const data = await response.json();
        setLigandSmiles(data.smiles);
      }
    } catch {
      // Fallback: just use content if it looks like SMILES
      if (!content.includes('\n') && content.length < 500) {
        setLigandSmiles(content.trim());
      }
    }
  };

  const calculateDockingBox = async (pdbContent: string) => {
    try {
      const response = await fetch('http://localhost:8000/api/v1/preparation/docking_box', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${localStorage.getItem('access_token')}`,
        },
        body: JSON.stringify({ pdb_content: pdbContent }),
      });
      
      if (response.ok) {
        const data = await response.json();
        setBoxCenter({
          x: data.center[0],
          y: data.center[1],
          z: data.center[2]
        });
        setBoxSize({
          x: data.size[0],
          y: data.size[1],
          z: data.size[2]
        });
      }
    } catch {
      // Ignore errors, keep default values
    }
  };

  const handleDocking = async () => {
    if (!proteinPdb.trim() || !ligandSmiles.trim()) return;

    setIsLoading(true);
    setError('');
    setPoses([]);

    try {
      const response = await fetch('http://localhost:8000/api/v1/docking/run', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${localStorage.getItem('access_token')}`,
        },
        body: JSON.stringify({
          protein_pdb: proteinPdb,
          ligand_smiles: ligandSmiles,
          box_center: [boxCenter.x, boxCenter.y, boxCenter.z],
          box_size: [boxSize.x, boxSize.y, boxSize.z],
        }),
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Docking failed');
      }

      const data = await response.json();
      setPoses(data.poses || []);
    } catch (err: any) {
      setError(err.message || 'Docking failed');
    } finally {
      setIsLoading(false);
    }
  };

  const getScoreColor = (score: number) => {
    if (score <= -8) return 'text-green-600';
    if (score <= -6) return 'text-yellow-600';
    return 'text-red-600';
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
              <Database className="h-8 w-8 text-purple-600" />
              <div>
                <h1 className="text-2xl font-bold text-slate-900">Molecular Docking</h1>
                <p className="text-sm text-slate-600">Protein-Ligand Binding Pose Prediction</p>
              </div>
            </div>
          </div>
        </div>
      </header>

      <main className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
          {/* Input Panel */}
          <div className="space-y-6">
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
                    placeholder="e.g., 1HSG, 6LU7"
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

              {/* Upload PDB */}
              <div className="mb-4">
                <label className="block text-sm font-medium text-slate-700 mb-1">
                  Or Upload PDB File
                </label>
                <input
                  type="file"
                  ref={proteinFileRef}
                  onChange={handleProteinFileUpload}
                  accept=".pdb,.pdbqt"
                  className="hidden"
                />
                <button
                  onClick={() => proteinFileRef.current?.click()}
                  className="w-full btn btn-secondary flex items-center justify-center space-x-2"
                >
                  <Upload className="h-4 w-4" />
                  <span>Upload PDB/PDBQT</span>
                </button>
              </div>

              {proteinPdb && (
                <div className="text-sm text-green-600 mb-4">
                  ✓ Protein loaded ({proteinPdb.split('\n').filter(l => l.startsWith('ATOM')).length} atoms)
                </div>
              )}
            </div>

            <div className="card">
              <h2 className="text-xl font-semibold mb-4">Ligand</h2>

              {/* Fetch from PubChem */}
              <div className="mb-4">
                <label className="block text-sm font-medium text-slate-700 mb-1">
                  Fetch from PubChem (CID)
                </label>
                <div className="flex space-x-2">
                  <input
                    type="text"
                    className="input flex-1"
                    value={pubchemCid}
                    onChange={(e) => setPubchemCid(e.target.value)}
                    placeholder="e.g., 2244 (Aspirin)"
                  />
                  <button
                    onClick={fetchFromPubChem}
                    className="btn btn-secondary flex items-center space-x-2"
                    disabled={isFetching || !pubchemCid.trim()}
                  >
                    <Search className="h-4 w-4" />
                    <span>Fetch</span>
                  </button>
                </div>
              </div>

              {/* Upload ligand file */}
              <div className="mb-4">
                <label className="block text-sm font-medium text-slate-700 mb-1">
                  Or Upload Ligand File
                </label>
                <input
                  type="file"
                  ref={ligandFileRef}
                  onChange={handleLigandFileUpload}
                  accept=".sdf,.mol,.mol2"
                  className="hidden"
                />
                <button
                  onClick={() => ligandFileRef.current?.click()}
                  className="w-full btn btn-secondary flex items-center justify-center space-x-2"
                >
                  <Upload className="h-4 w-4" />
                  <span>Upload SDF/MOL/MOL2</span>
                </button>
              </div>

              <div>
                <label className="block text-sm font-medium text-slate-700 mb-1">
                  Or Enter SMILES
                </label>
                <input
                  type="text"
                  className="input font-mono"
                  value={ligandSmiles}
                  onChange={(e) => setLigandSmiles(e.target.value)}
                  placeholder="e.g., CC(=O)Oc1ccccc1C(=O)O"
                />
              </div>
            </div>

            <div className="card">
              <h2 className="text-xl font-semibold mb-4">Docking Box</h2>
              
              <div className="grid grid-cols-3 gap-4 mb-4">
                <div>
                  <label className="text-xs text-slate-500">Center X</label>
                  <input
                    type="number"
                    className="input"
                    value={boxCenter.x}
                    onChange={(e) => setBoxCenter({ ...boxCenter, x: parseFloat(e.target.value) || 0 })}
                  />
                </div>
                <div>
                  <label className="text-xs text-slate-500">Center Y</label>
                  <input
                    type="number"
                    className="input"
                    value={boxCenter.y}
                    onChange={(e) => setBoxCenter({ ...boxCenter, y: parseFloat(e.target.value) || 0 })}
                  />
                </div>
                <div>
                  <label className="text-xs text-slate-500">Center Z</label>
                  <input
                    type="number"
                    className="input"
                    value={boxCenter.z}
                    onChange={(e) => setBoxCenter({ ...boxCenter, z: parseFloat(e.target.value) || 0 })}
                  />
                </div>
              </div>

              <div className="grid grid-cols-3 gap-4">
                <div>
                  <label className="text-xs text-slate-500">Size X (Å)</label>
                  <input
                    type="number"
                    className="input"
                    value={boxSize.x}
                    onChange={(e) => setBoxSize({ ...boxSize, x: parseFloat(e.target.value) || 20 })}
                  />
                </div>
                <div>
                  <label className="text-xs text-slate-500">Size Y (Å)</label>
                  <input
                    type="number"
                    className="input"
                    value={boxSize.y}
                    onChange={(e) => setBoxSize({ ...boxSize, y: parseFloat(e.target.value) || 20 })}
                  />
                </div>
                <div>
                  <label className="text-xs text-slate-500">Size Z (Å)</label>
                  <input
                    type="number"
                    className="input"
                    value={boxSize.z}
                    onChange={(e) => setBoxSize({ ...boxSize, z: parseFloat(e.target.value) || 20 })}
                  />
                </div>
              </div>

              <button
                onClick={handleDocking}
                className="w-full btn btn-primary mt-6"
                disabled={isLoading || !proteinPdb.trim() || !ligandSmiles.trim()}
              >
                {isLoading ? (
                  <span className="flex items-center justify-center space-x-2">
                    <Loader className="h-4 w-4 animate-spin" />
                    <span>Running Docking...</span>
                  </span>
                ) : (
                  'Run Docking'
                )}
              </button>

              {error && (
                <div className="mt-4 bg-red-50 border border-red-200 text-red-700 px-4 py-3 rounded-lg text-sm">
                  {error}
                </div>
              )}
            </div>
          </div>

          {/* Results Panel */}
          <div className="card">
            <h2 className="text-xl font-semibold mb-4">Docking Results</h2>

            {poses.length === 0 && !error && (
              <div className="text-center py-12 text-slate-500">
                <Database className="h-16 w-16 mx-auto mb-4 text-slate-300" />
                <p>Load a protein and ligand, then run docking</p>
              </div>
            )}

            {poses.length > 0 && (
              <div className="space-y-4">
                <div className="bg-purple-50 border border-purple-200 rounded-lg p-4">
                  <p className="font-medium text-purple-900">
                    {poses.length} pose(s) generated
                  </p>
                  <p className="text-sm text-purple-700 mt-1">
                    Best score: {poses[0]?.score?.toFixed(2)} kcal/mol
                  </p>
                </div>

                <div className="space-y-2">
                  {poses.map((pose, idx) => (
                    <div
                      key={idx}
                      className="flex items-center justify-between p-4 bg-slate-50 rounded-lg"
                    >
                      <div>
                        <p className="font-medium">Pose {idx + 1}</p>
                        <p className="text-sm text-slate-600">
                          RMSD: {pose.rmsd?.toFixed(2)} Å
                        </p>
                      </div>
                      <div className="text-right">
                        <p className={`text-lg font-bold ${getScoreColor(pose.score)}`}>
                          {pose.score?.toFixed(2)} kcal/mol
                        </p>
                        <p className="text-xs text-slate-500">
                          LE: {pose.ligand_efficiency?.toFixed(2)}
                        </p>
                      </div>
                    </div>
                  ))}
                </div>

                <div className="bg-yellow-50 border border-yellow-200 text-yellow-800 px-4 py-3 rounded-lg text-xs">
                  <p className="font-medium">Note:</p>
                  <p>Docking scores are predicted binding affinity estimates, not experimental values.</p>
                </div>
              </div>
            )}
          </div>
        </div>
      </main>
    </div>
  );
};
