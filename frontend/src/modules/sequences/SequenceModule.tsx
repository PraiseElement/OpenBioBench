/**
 * Sequence Alignment Module.
 * Multiple sequence alignment with database fetch and conservation visualization.
 */
import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import { AlignLeft, ArrowLeft, Loader, Download, Search } from 'lucide-react';

export const SequenceModule: React.FC = () => {
  const [fastaInput, setFastaInput] = useState('');
  const [alignmentResult, setAlignmentResult] = useState<any>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState('');
  
  // Database search state
  const [searchQuery, setSearchQuery] = useState('');
  const [searchResults, setSearchResults] = useState<any[]>([]);
  const [isSearching, setIsSearching] = useState(false);
  const [showSearchModal, setShowSearchModal] = useState(false);

  const exampleFasta = `>Sequence1_Human
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH
>Sequence2_Mouse
MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSH
>Sequence3_Chicken
MVLSAADKNNVKGIFTKIAGHAEEYGAETLERMFTTYPPTKTYFPHFDLSH`;

  const searchUniProt = async () => {
    if (!searchQuery.trim()) return;
    setIsSearching(true);
    setSearchResults([]);
    
    try {
      const response = await fetch(
        `http://localhost:8000/api/v1/database/uniprot/search/?query=${encodeURIComponent(searchQuery)}&limit=10`,
        { headers: { 'Authorization': `Bearer ${localStorage.getItem('access_token')}` } }
      );
      
      if (!response.ok) throw new Error('Search failed');
      
      const data = await response.json();
      setSearchResults(data);
    } catch (err) {
      setError('Search failed');
    } finally {
      setIsSearching(false);
    }
  };

  const fetchSequence = async (accession: string) => {
    try {
      const response = await fetch(
        `http://localhost:8000/api/v1/database/uniprot/${accession}`,
        { headers: { 'Authorization': `Bearer ${localStorage.getItem('access_token')}` } }
      );
      
      if (!response.ok) throw new Error('Fetch failed');
      
      const data = await response.json();
      const newFasta = `>${accession}\n${data.sequence}`;
      setFastaInput(prev => prev ? `${prev}\n\n${newFasta}` : newFasta);
      setShowSearchModal(false);
    } catch (err) {
      setError('Failed to fetch sequence');
    }
  };

  const handleAlign = async () => {
    if (!fastaInput.trim()) return;

    setIsLoading(true);
    setError('');
    setAlignmentResult(null);

    try {
      const response = await fetch('http://localhost:8000/api/v1/alignments/align', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
          'Authorization': `Bearer ${localStorage.getItem('access_token')}`,
        },
        body: JSON.stringify({ fasta_content: fastaInput }),
      });

      if (!response.ok) {
        const errorData = await response.json();
        throw new Error(errorData.detail || 'Alignment failed');
      }

      const data = await response.json();
      setAlignmentResult(data);
    } catch (err: any) {
      setError(err.message || 'Alignment failed');
    } finally {
      setIsLoading(false);
    }
  };

  const getConservationColor = (conservation: number) => {
    if (conservation >= 0.9) return 'bg-blue-600 text-white';
    if (conservation >= 0.7) return 'bg-blue-400 text-white';
    if (conservation >= 0.5) return 'bg-blue-200';
    return 'bg-slate-100';
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
              <AlignLeft className="h-8 w-8 text-cyan-600" />
              <div>
                <h1 className="text-2xl font-bold text-slate-900 dark:text-slate-100">Sequence Alignment</h1>
                <p className="text-sm text-slate-600 dark:text-slate-400">Multiple Sequence Alignment & Conservation</p>
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
              <h2 className="text-xl font-semibold mb-4">Input Sequences</h2>

              {/* UniProt Search Button */}
              <button
                onClick={() => setShowSearchModal(true)}
                className="w-full btn btn-secondary mb-4 flex items-center justify-center space-x-2"
              >
                <Search className="h-4 w-4" />
                <span>Search UniProt Database</span>
              </button>

              <div className="space-y-4">
                <div>
                  <label className="block text-sm font-medium text-slate-700 mb-1">
                    FASTA Sequences (2+ sequences required)
                  </label>
                  <textarea
                    className="input font-mono text-sm"
                    rows={12}
                    value={fastaInput}
                    onChange={(e) => setFastaInput(e.target.value)}
                    placeholder={exampleFasta}
                  />
                </div>

                <button
                  onClick={() => setFastaInput(exampleFasta)}
                  className="text-sm text-cyan-600 hover:text-cyan-700"
                >
                  Load example sequences
                </button>

                <button
                  onClick={handleAlign}
                  className="w-full btn btn-primary"
                  disabled={isLoading || !fastaInput.trim()}
                >
                  {isLoading ? (
                    <span className="flex items-center justify-center space-x-2">
                      <Loader className="h-4 w-4 animate-spin" />
                      <span>Aligning...</span>
                    </span>
                  ) : (
                    'Run Alignment'
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
            <h2 className="text-xl font-semibold mb-4">Alignment Results</h2>

            {!alignmentResult && !error && (
              <div className="text-center py-12 text-slate-500">
                <AlignLeft className="h-16 w-16 mx-auto mb-4 text-slate-300" />
                <p>Enter FASTA sequences and run alignment to see results</p>
              </div>
            )}

            {alignmentResult && (
              <div className="space-y-4">
                <div className="bg-cyan-50 border border-cyan-200 rounded-lg p-4">
                  <h3 className="font-semibold text-cyan-900 mb-2">Alignment Statistics</h3>
                  <div className="grid grid-cols-2 gap-2 text-sm">
                    <div>Sequences: {alignmentResult.num_sequences}</div>
                    <div>Length: {alignmentResult.alignment_length}</div>
                    <div>Identity: {(alignmentResult.percent_identity * 100).toFixed(1)}%</div>
                    <div>Gaps: {(alignmentResult.gap_percentage * 100).toFixed(1)}%</div>
                  </div>
                </div>

                <div>
                  <h3 className="font-semibold mb-2">Aligned Sequences</h3>
                  <div className="overflow-x-auto">
                    <div className="font-mono text-xs space-y-1 bg-slate-50 p-4 rounded-lg">
                      {alignmentResult.aligned_sequences?.map((seq: any, idx: number) => (
                        <div key={idx} className="flex">
                          <span className="w-24 truncate text-slate-600 mr-2">{seq.id}:</span>
                          <span className="flex">
                            {seq.sequence.split('').map((char: string, i: number) => (
                              <span
                                key={i}
                                className={`w-4 text-center ${
                                  alignmentResult.conservation?.[i]
                                    ? getConservationColor(alignmentResult.conservation[i])
                                    : ''
                                }`}
                              >
                                {char}
                              </span>
                            ))}
                          </span>
                        </div>
                      ))}
                    </div>
                  </div>
                </div>

                <div className="flex space-x-2">
                  <button className="btn btn-secondary flex items-center space-x-2">
                    <Download className="h-4 w-4" />
                    <span>Export FASTA</span>
                  </button>
                </div>
              </div>
            )}
          </div>
        </div>
      </main>

      {/* UniProt Search Modal */}
      {showSearchModal && (
        <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center z-50">
          <div className="bg-white rounded-xl p-6 w-full max-w-2xl max-h-[80vh] overflow-y-auto">
            <h3 className="text-xl font-semibold mb-4">Search UniProt Database</h3>
            
            <div className="flex space-x-2 mb-4">
              <input
                type="text"
                className="input flex-1"
                value={searchQuery}
                onChange={(e) => setSearchQuery(e.target.value)}
                placeholder="Search by gene name, protein name, organism..."
                onKeyDown={(e) => e.key === 'Enter' && searchUniProt()}
              />
              <button
                onClick={searchUniProt}
                className="btn btn-primary"
                disabled={isSearching}
              >
                {isSearching ? <Loader className="h-4 w-4 animate-spin" /> : 'Search'}
              </button>
            </div>

            <div className="text-sm text-slate-500 mb-4">
              Examples: EGFR, p53, hemoglobin human, kinase
            </div>

            {searchResults.length > 0 && (
              <div className="space-y-2">
                {searchResults.map((result) => (
                  <div
                    key={result.accession}
                    className="p-3 border rounded-lg hover:bg-slate-50 cursor-pointer"
                    onClick={() => fetchSequence(result.accession)}
                  >
                    <div className="font-medium">{result.accession} - {result.protein_name}</div>
                    <div className="text-sm text-slate-600">
                      {result.organism} | {result.length} aa
                    </div>
                  </div>
                ))}
              </div>
            )}

            <div className="mt-4 flex justify-end">
              <button
                onClick={() => setShowSearchModal(false)}
                className="btn btn-secondary"
              >
                Close
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};
