/**
 * ADMET Prediction Module.
 */
import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import { apiClient } from '../../services/api';
import { Activity, ArrowLeft, Loader, AlertTriangle, CheckCircle } from 'lucide-react';
import type { ADMETResponse } from '../../types';

export const ADMETModule: React.FC = () => {
  const [smiles, setSmiles] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [result, setResult] = useState<ADMETResponse | null>(null);
  const [error, setError] = useState('');

  const handlePredict = async () => {
    if (!smiles.trim()) return;

    setIsLoading(true);
    setError('');
    setResult(null);

    try {
      const response = await apiClient.predictADMET({ smiles });
      setResult(response);
    } catch (err: any) {
      setError(err.response?.data?.detail || 'Prediction failed');
    } finally {
      setIsLoading(false);
    }
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
              <Activity className="h-8 w-8 text-green-600" />
              <div>
                <h1 className="text-2xl font-bold text-slate-900 dark:text-slate-100">ADMET Prediction</h1>
                <p className="text-sm text-slate-600 dark:text-slate-400">Pharmacokinetic & Toxicity Properties</p>
              </div>
            </div>
          </div>
        </div>
      </header>

      <main className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
          {/* Input Panel */}
          <div className="lg:col-span-1">
            <div className="card sticky top-8">
              <h2 className="text-xl font-semibold mb-4">Molecule Input</h2>

              <div className="space-y-4">
                <div>
                  <label className="block text-sm font-medium text-slate-700 mb-1">
                    SMILES String *
                  </label>
                  <textarea
                    className="input font-mono text-sm"
                    rows={4}
                    value={smiles}
                    onChange={(e) => setSmiles(e.target.value)}
                    placeholder="CC(=O)Oc1ccccc1C(=O)O"
                  />
                </div>

                <button
                  onClick={handlePredict}
                  className="w-full btn btn-primary"
                  disabled={isLoading || !smiles.trim()}
                >
                  {isLoading ? (
                    <span className="flex items-center justify-center space-x-2">
                      <Loader className="h-4 w-4 animate-spin" />
                      <span>Predicting...</span>
                    </span>
                  ) : (
                    'Predict Properties'
                  )}
                </button>
              </div>

              {error && (
                <div className="mt-4 bg-red-50 border border-red-200 text-red-700 px-4 py-3 rounded-lg text-sm">
                  {error}
                </div>
              )}
            </div>
          </div>

          {/* Results Panel */}
          <div className="lg:col-span-2 space-y-6">
            {!result && !error && (
              <div className="card text-center py-12 text-slate-500">
                <Activity className="h-16 w-16 mx-auto mb-4 text-slate-300" />
                <p>Enter a SMILES string to predict ADMET properties</p>
              </div>
            )}

            {result && (
              <>
                {/* Physicochemical Properties */}
                <div className="card">
                  <h3 className="text-lg font-semibold mb-4">Physicochemical Properties</h3>
                  <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                    <div className="bg-slate-50 p-3 rounded">
                      <p className="text-xs text-slate-600 mb-1">Molecular Weight</p>
                      <p className="text-lg font-semibold">{result.properties.molecular_weight.toFixed(2)}</p>
                      <p className="text-xs text-slate-500">g/mol</p>
                    </div>
                    <div className="bg-slate-50 p-3 rounded">
                      <p className="text-xs text-slate-600 mb-1">LogP</p>
                      <p className="text-lg font-semibold">{result.properties.logp.toFixed(2)}</p>
                    </div>
                    <div className="bg-slate-50 p-3 rounded">
                      <p className="text-xs text-slate-600 mb-1">TPSA</p>
                      <p className="text-lg font-semibold">{result.properties.tpsa.toFixed(1)}</p>
                      <p className="text-xs text-slate-500">Ų</p>
                    </div>
                    <div className="bg-slate-50 p-3 rounded">
                      <p className="text-xs text-slate-600 mb-1">Rotatable Bonds</p>
                      <p className="text-lg font-semibold">{result.properties.rotatable_bonds}</p>
                    </div>
                  </div>

                  <div className="mt-4 flex items-center space-x-2">
                    {result.properties.lipinski_violations === 0 ? (
                      <>
                        <CheckCircle className="h-5 w-5 text-green-600" />
                        <span className="text-sm text-green-700 font-medium">
                          Lipinski Rule of Five: Passed
                        </span>
                      </>
                    ) : (
                      <>
                        <AlertTriangle className="h-5 w-5 text-orange-600" />
                        <span className="text-sm text-orange-700 font-medium">
                          {result.properties.lipinski_violations} Lipinski violation(s)
                        </span>
                      </>
                    )}
                  </div>
                </div>

                {/* Pharmacokinetics */}
                <div className="card">
                  <h3 className="text-lg font-semibold mb-4">Pharmacokinetics</h3>
                  <div className="space-y-3">
                    <div className="bg-slate-50 p-4 rounded">
                      <div className="flex justify-between items-start mb-2">
                        <div>
                          <p className="font-medium text-slate-900">Caco-2 Permeability</p>
                          <p className="text-sm text-slate-600 mt-1">
                            {result.pharmacokinetics.caco2_permeability.interpretation}
                          </p>
                        </div>
                        <span className="text-xs px-2 py-1 bg-slate-200 rounded">
                          {result.pharmacokinetics.caco2_permeability.confidence}
                        </span>
                      </div>
                      <p className="text-sm text-slate-500">
                        {result.pharmacokinetics.caco2_permeability.value} {result.pharmacokinetics.caco2_permeability.unit}
                      </p>
                    </div>

                    <div className="bg-slate-50 p-4 rounded">
                      <div className="flex justify-between items-start mb-2">
                        <div>
                          <p className="font-medium text-slate-900">BBB Penetration</p>
                          <p className="text-sm text-slate-600 mt-1">
                            {result.pharmacokinetics.bbb_penetration.interpretation}
                          </p>
                        </div>
                        <span className="text-xs px-2 py-1 bg-slate-200 rounded">
                          {result.pharmacokinetics.bbb_penetration.confidence}
                        </span>
                      </div>
                      <p className="text-sm text-slate-500">
                        Score: {result.pharmacokinetics.bbb_penetration.value}
                      </p>
                    </div>
                  </div>
                </div>

                {/* Toxicity */}
                <div className="card">
                  <h3 className="text-lg font-semibold mb-4">Toxicity</h3>
                  <div className="space-y-3">
                    <div className="bg-slate-50 p-4 rounded">
                      <div className="flex justify-between items-start mb-2">
                        <div>
                          <p className="font-medium text-slate-900">hERG Cardiotoxicity</p>
                          <p className="text-sm text-slate-600 mt-1">
                            Risk: {result.toxicity.herg_ic50.interpretation || 'Unknown'}
                          </p>
                        </div>
                        <span className="text-xs px-2 py-1 bg-slate-200 rounded">
                          {result.toxicity.herg_ic50.confidence}
                        </span>
                      </div>
                      <p className="text-sm text-slate-500">
                        IC50: {result.toxicity.herg_ic50.value} {result.toxicity.herg_ic50.unit}
                      </p>
                    </div>

                    <div className="bg-slate-50 p-4 rounded">
                      <div className="flex justify-between items-start">
                        <div>
                          <p className="font-medium text-slate-900">AMES Mutagenicity</p>
                          <p className="text-sm text-slate-600 mt-1">
                            {result.toxicity.ames_mutagenicity.prediction}
                          </p>
                        </div>
                        <span className="text-xs px-2 py-1 bg-slate-200 rounded">
                          Confidence: {(result.toxicity.ames_mutagenicity.confidence * 100).toFixed(0)}%
                        </span>
                      </div>
                    </div>
                  </div>
                </div>

                {/* Warnings */}
                {result.warnings.length > 0 && (
                  <div className="bg-yellow-50 border border-yellow-200 rounded-lg p-4">
                    <div className="flex items-start space-x-3">
                      <AlertTriangle className="h-5 w-5 text-yellow-600 flex-shrink-0 mt-0.5" />
                      <div className="flex-1">
                        <p className="font-medium text-yellow-900 mb-2">Important Notes:</p>
                        <ul className="text-sm text-yellow-800 space-y-1">
                          {result.warnings.map((warning, idx) => (
                            <li key={idx}>• {warning}</li>
                          ))}
                        </ul>
                      </div>
                    </div>
                  </div>
                )}
              </>
            )}
          </div>
        </div>
      </main>
    </div>
  );
};
