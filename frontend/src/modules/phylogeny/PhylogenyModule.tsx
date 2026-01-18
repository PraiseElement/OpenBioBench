/**
 * Phylogenetics Module.
 * Phylogenetic tree construction from alignments.
 */
import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import { GitBranch, ArrowLeft, Loader, Download } from 'lucide-react';

export const PhylogenyModule: React.FC = () => {
  const [alignmentInput, setAlignmentInput] = useState('');
  const [method, setMethod] = useState('nj');
  const [treeResult, setTreeResult] = useState<any>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState('');

  const exampleAlignment = `>Human
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH
>Mouse
MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSH
>Chicken
MVLSAADKNNVKGIFTKIAGHAEEYGAETLERMFTTYPPTKTYFPHFDLSH
>Zebrafish
MSLSDKDKAAVRALWSKIGKSADAIGNDALSRMIVVYPQTKTYFSHWADLS`;

  const handleBuildTree = async () => {
    if (!alignmentInput.trim()) return;

    setIsLoading(true);
    setError('');
    setTreeResult(null);

    try {
      const response = await fetch('http://localhost:8000/api/v1/phylogeny/infer', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${localStorage.getItem('access_token')}`,
        },
        body: JSON.stringify({
          alignment_fasta: alignmentInput,
          method: method,
        }),
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Tree construction failed');
      }

      const data = await response.json();
      setTreeResult(data);
    } catch (err: any) {
      setError(err.message || 'Tree construction failed');
    } finally {
      setIsLoading(false);
    }
  };

  const renderAsciiTree = (newick: string) => {
    // Simple ASCII representation - in production would use a proper tree viz library
    return (
      <pre className="font-mono text-sm bg-slate-900 text-green-400 p-4 rounded-lg overflow-x-auto">
        {newick || 'No tree data available'}
      </pre>
    );
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
              <GitBranch className="h-8 w-8 text-teal-600" />
              <div>
                <h1 className="text-2xl font-bold text-slate-900 dark:text-slate-100">Phylogenetic Analysis</h1>
                <p className="text-sm text-slate-600 dark:text-slate-400">Evolutionary Tree Construction</p>
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
              <h2 className="text-xl font-semibold mb-4">Input Alignment</h2>

              <div className="space-y-4">
                <div>
                  <label className="block text-sm font-medium text-slate-700 mb-1">
                    Aligned Sequences (FASTA format)
                  </label>
                  <textarea
                    className="input font-mono text-sm"
                    rows={10}
                    value={alignmentInput}
                    onChange={(e) => setAlignmentInput(e.target.value)}
                    placeholder={exampleAlignment}
                  />
                </div>

                <button
                  onClick={() => setAlignmentInput(exampleAlignment)}
                  className="text-sm text-teal-600 hover:text-teal-700"
                >
                  Load example alignment
                </button>

                <div>
                  <label className="block text-sm font-medium text-slate-700 mb-1">
                    Tree Construction Method
                  </label>
                  <select
                    className="input"
                    value={method}
                    onChange={(e) => setMethod(e.target.value)}
                  >
                    <option value="nj">Neighbor-Joining (NJ)</option>
                    <option value="upgma">UPGMA</option>
                  </select>
                </div>

                <button
                  onClick={handleBuildTree}
                  className="w-full btn btn-primary"
                  disabled={isLoading || !alignmentInput.trim()}
                >
                  {isLoading ? (
                    <span className="flex items-center justify-center space-x-2">
                      <Loader className="h-4 w-4 animate-spin" />
                      <span>Building Tree...</span>
                    </span>
                  ) : (
                    'Build Phylogenetic Tree'
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
          <div className="card">
            <h2 className="text-xl font-semibold mb-4">Phylogenetic Tree</h2>

            {!treeResult && !error && (
              <div className="text-center py-12 text-slate-500">
                <GitBranch className="h-16 w-16 mx-auto mb-4 text-slate-300" />
                <p>Enter aligned sequences to build a phylogenetic tree</p>
              </div>
            )}

            {treeResult && (
              <div className="space-y-4">
                <div className="bg-teal-50 border border-teal-200 rounded-lg p-4">
                  <h3 className="font-semibold text-teal-900 mb-2">Tree Statistics</h3>
                  <div className="grid grid-cols-2 gap-2 text-sm">
                    <div>Taxa: {treeResult.num_taxa}</div>
                    <div>Method: {treeResult.method?.toUpperCase()}</div>
                    <div>Total branch length: {treeResult.total_branch_length?.toFixed(4)}</div>
                  </div>
                </div>

                <div>
                  <h3 className="font-semibold mb-2">Newick Format</h3>
                  {renderAsciiTree(treeResult.newick)}
                </div>

                {treeResult.ascii_tree && (
                  <div>
                    <h3 className="font-semibold mb-2">Tree Visualization</h3>
                    <pre className="font-mono text-xs bg-slate-50 p-4 rounded-lg overflow-x-auto whitespace-pre">
                      {treeResult.ascii_tree}
                    </pre>
                  </div>
                )}

                <div className="flex space-x-2">
                  <button className="btn btn-secondary flex items-center space-x-2">
                    <Download className="h-4 w-4" />
                    <span>Export Newick</span>
                  </button>
                </div>
              </div>
            )}
          </div>
        </div>
      </main>
    </div>
  );
};
