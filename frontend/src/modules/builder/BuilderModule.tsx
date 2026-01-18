/**
 * Ligand Builder Module - 2D/3D molecular design.
 */
import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import { apiClient } from '../../services/api';
import { FlaskConical, ArrowLeft, Loader } from 'lucide-react';
import type { Generate3DResponse } from '../../types';

export const BuilderModule: React.FC = () => {
  const [smiles, setSmiles] = useState('');
  const [name, setName] = useState('');
  const [isValidating, setIsValidating] = useState(false);
  const [isGenerating, setIsGenerating] = useState(false);
  const [validationResult, setValidationResult] = useState<any>(null);
  const [result, setResult] = useState<Generate3DResponse | null>(null);
  const [error, setError] = useState('');

  // Hard-coded project ID for demo (in production, user would select project)
  const demoProjectId = '00000000-0000-0000-0000-000000000001';

  const handleValidate = async () => {
    if (!smiles.trim()) return;

    setIsValidating(true);
    setError('');
    setValidationResult(null);

    try {
      const result = await apiClient.validateSMILES(smiles);
      setValidationResult(result);
      if (result.valid) {
        setSmiles(result.canonical_smiles);
      }
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Validation failed');
    } finally {
      setIsValidating(false);
    }
  };

  const handleGenerate3D = async () => {
    setIsGenerating(true);
    setError('');
    setResult(null);

    try {
      const response = await apiClient.generate3D({
        smiles,
        project_id: demoProjectId,
        num_conformers: 10,
        minimize: true,
        force_field: 'mmff94',
        name: name || undefined,
      });
      setResult(response);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Generation failed');
    } finally {
      setIsGenerating(false);
    }
  };

  return (
    <div className="min-h-screen bg-slate-50">
      {/* Header */}
      <header className="bg-white border-b">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-4">
          <div className="flex items-center justify-between">
            <div className="flex items-center space-x-4">
              <Link to="/" className="text-slate-600 hover:text-slate-900">
                <ArrowLeft className="h-6 w-6" />
              </Link>
              <div className="flex items-center space-x-3">
                <FlaskConical className="h-8 w-8 text-primary-600" />
                <div>
                  <h1 className="text-2xl font-bold text-slate-900">Ligand Builder</h1>
                  <p className="text-sm text-slate-600">2D/3D Molecular Design</p>
                </div>
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
              <h2 className="text-xl font-semibold mb-4">Molecule Input</h2>

              <div className="space-y-4">
                <div>
                  <label className="block text-sm font-medium text-slate-700 mb-1">
                    SMILES String *
                  </label>
                  <textarea
                    className="input font-mono text-sm"
                    rows={3}
                    value={smiles}
                    onChange={(e) => setSmiles(e.target.value)}
                    placeholder="CC(=O)Oc1ccccc1C(=O)O"
                  />
                  <p className="mt-1 text-xs text-slate-500">
                    Example: Aspirin = CC(=O)Oc1ccccc1C(=O)O
                  </p>
                </div>

                <div>
                  <label className="block text-sm font-medium text-slate-700 mb-1">
                    Molecule Name (Optional)
                  </label>
                  <input
                    type="text"
                    className="input"
                    value={name}
                    onChange={(e) => setName(e.target.value)}
                    placeholder="e.g., Aspirin"
                  />
                </div>

                <button
                  onClick={handleValidate}
                  className="w-full btn btn-secondary"
                  disabled={isValidating || !smiles.trim()}
                >
                  {isValidating ? 'Validating...' : 'Validate SMILES'}
                </button>

                {validationResult && validationResult.valid && (
                  <button
                    onClick={handleGenerate3D}
                    className="w-full btn btn-primary"
                    disabled={isGenerating}
                  >
                    {isGenerating ? (
                      <span className="flex items-center justify-center space-x-2">
                        <Loader className="h-4 w-4 animate-spin" />
                        <span>Generating 3D...</span>
                      </span>
                    ) : (
                      'Generate 3D Conformers'
                    )}
                  </button>
                )}
              </div>

              {error && (
                <div className="mt-4 bg-red-50 border border-red-200 text-red-700 px-4 py-3 rounded-lg text-sm">
                  {error}
                </div>
              )}

              {validationResult && validationResult.valid && (
                <div className="mt-4 bg-green-50 border border-green-200 text-green-700 px-4 py-3 rounded-lg">
                  <p className="font-medium text-sm mb-2">✓ Valid SMILES</p>
                  <div className="text-xs space-y-1">
                    <p>Formula: {validationResult.molecular_formula}</p>
                    <p>MW: {validationResult.molecular_weight?.toFixed(2)} g/mol</p>
                  </div>
                </div>
              )}
            </div>
          </div>

          {/* Results Panel */}
          <div className="card">
            <h2 className="text-xl font-semibold mb-4">Results</h2>

            {!result && (
              <div className="text-center py-12 text-slate-500">
                <p>Enter a SMILES string and generate 3D conformers to see results</p>
              </div>
            )}

            {result && (
              <div className="space-y-4">
                <div className="bg-primary-50 border border-primary-200 rounded-lg p-4">
                  <h3 className="font-semibold text-primary-900 mb-2">Generated</h3>
                  <p className="text-sm text-primary-700 mb-2">
                    {result.conformers.length} conformers
                  </p>
                  <p className="text-xs font-mono text-primary-600">
                    {result.canonical_smiles}
                  </p>
                </div>

                <div>
                  <h3 className="font-semibold mb-2">Molecular Properties</h3>
                  <div className="grid grid-cols-2 gap-2 text-sm">
                    <div className="bg-slate-50 p-2 rounded">
                      <p className="text-slate-600">MW</p>
                      <p className="font-semibold">{result.properties.molecular_weight.toFixed(2)}</p>
                    </div>
                    <div className="bg-slate-50 p-2 rounded">
                      <p className="text-slate-600">LogP</p>
                      <p className="font-semibold">{result.properties.logp.toFixed(2)}</p>
                    </div>
                    <div className="bg-slate-50 p-2 rounded">
                      <p className="text-slate-600">HBD</p>
                      <p className="font-semibold">{result.properties.hbd}</p>
                    </div>
                    <div className="bg-slate-50 p-2 rounded">
                      <p className="text-slate-600">HBA</p>
                      <p className="font-semibold">{result.properties.hba}</p>
                    </div>
                  </div>
                </div>

                <div>
                  <h3 className="font-semibold mb-2">Conformers</h3>
                  <div className="max-h-64 overflow-y-auto space-y-2">
                    {result.conformers.map((conf) => (
                      <div
                        key={conf.conformer_id}
                        className="bg-slate-50 p-3 rounded text-sm"
                      >
                        <div className="flex justify-between items-center">
                          <span className="font-medium">Conformer {conf.conformer_id}</span>
                          <span className="text-xs text-slate-600">
                            {conf.energy.toFixed(2)} kcal/mol
                          </span>
                        </div>
                        {conf.rmsd !== null && conf.rmsd !== undefined && (
                          <p className="text-xs text-slate-500 mt-1">
                            RMSD: {conf.rmsd.toFixed(3)} Å
                          </p>
                        )}
                      </div>
                    ))}
                  </div>
                </div>

                <div className="bg-yellow-50 border border-yellow-200 text-yellow-800 px-4 py-3 rounded-lg text-xs">
                  <p className="font-medium">Note:</p>
                  <p>These are computational predictions for research purposes only.</p>
                </div>
              </div>
            )}
          </div>
        </div>
      </main>
    </div>
  );
};
