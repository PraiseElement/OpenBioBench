/**
 * Dashboard - Main landing page after login.
 * Features scientific aesthetic with molecular pattern background.
 */
import React from 'react';
import { Link } from 'react-router-dom';
import { useAuthStore } from '../../store/useAuthStore';
import { useThemeStore } from '../../store/useThemeStore';
import { 
  FlaskConical, 
  Beaker, 
  Activity, 
  Database, 
  LogOut, 
  AlignLeft, 
  Target, 
  GitBranch, 
  Wrench,
  Settings,
  Sun,
  Moon
} from 'lucide-react';

export const Dashboard: React.FC = () => {
  const { user, logout } = useAuthStore();
  const { theme, toggleTheme } = useThemeStore();

  const modules = [
    {
      name: 'Ligand Builder',
      description: '2D/3D molecular design with conformer generation',
      icon: Beaker,
      path: '/builder',
      color: 'bg-blue-500',
      gradient: 'from-blue-500 to-blue-600',
    },
    {
      name: 'ADMET Prediction',
      description: 'Predict pharmacokinetic and toxicity properties',
      icon: Activity,
      path: '/admet',
      color: 'bg-green-500',
      gradient: 'from-green-500 to-emerald-600',
    },
    {
      name: 'Sequence Alignment',
      description: 'Multiple sequence alignment with conservation analysis',
      icon: AlignLeft,
      path: '/sequences',
      color: 'bg-cyan-500',
      gradient: 'from-cyan-500 to-teal-600',
    },
    {
      name: 'Pocket Detection',
      description: 'Binding site identification and druggability scoring',
      icon: Target,
      path: '/pockets',
      color: 'bg-orange-500',
      gradient: 'from-orange-500 to-amber-600',
    },
    {
      name: 'Molecular Docking',
      description: 'Protein-ligand docking with pose prediction',
      icon: Database,
      path: '/docking',
      color: 'bg-purple-500',
      gradient: 'from-purple-500 to-violet-600',
    },
    {
      name: 'Phylogenetics',
      description: 'Phylogenetic tree construction from alignments',
      icon: GitBranch,
      path: '/phylogeny',
      color: 'bg-teal-500',
      gradient: 'from-teal-500 to-cyan-600',
    },
    {
      name: 'Structure Preparation',
      description: 'Clean & prepare proteins and ligands for docking',
      icon: Wrench,
      path: '/preparation',
      color: 'bg-amber-500',
      gradient: 'from-amber-500 to-yellow-600',
    },
  ];

  return (
    <div className={`min-h-screen bg-molecular ${theme === 'dark' ? 'bg-slate-900' : 'bg-slate-50'}`}>
      {/* Header */}
      <header className={`${theme === 'dark' ? 'bg-slate-800/80 border-slate-700' : 'bg-white/80'} backdrop-blur-md border-b sticky top-0 z-50`}>
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-4">
          <div className="flex items-center justify-between">
            <div className="flex items-center space-x-3">
              <div className="relative">
                <FlaskConical className="h-10 w-10 text-blue-500" />
                <div className="absolute -top-1 -right-1 w-3 h-3 bg-green-400 rounded-full pulse-ring" />
              </div>
              <div>
                <h1 className={`text-2xl font-bold ${theme === 'dark' ? 'text-white' : 'text-slate-900'}`}>
                  OpenBioBench
                </h1>
                <p className={`text-sm ${theme === 'dark' ? 'text-slate-400' : 'text-slate-600'}`}>
                  Computational Biology Platform
                </p>
              </div>
            </div>
            <div className="flex items-center space-x-4">
              {/* Theme Toggle */}
              <button
                onClick={toggleTheme}
                className={`p-2 rounded-lg transition-colors ${
                  theme === 'dark' 
                    ? 'bg-slate-700 text-yellow-400 hover:bg-slate-600' 
                    : 'bg-slate-100 text-slate-600 hover:bg-slate-200'
                }`}
                title="Toggle theme"
              >
                {theme === 'dark' ? <Sun className="h-5 w-5" /> : <Moon className="h-5 w-5" />}
              </button>
              
              {/* Settings */}
              <Link
                to="/settings"
                className={`p-2 rounded-lg transition-colors ${
                  theme === 'dark' 
                    ? 'bg-slate-700 text-slate-300 hover:bg-slate-600' 
                    : 'bg-slate-100 text-slate-600 hover:bg-slate-200'
                }`}
                title="Settings"
              >
                <Settings className="h-5 w-5" />
              </Link>

              <span className={`${theme === 'dark' ? 'text-slate-300' : 'text-slate-700'}`}>
                {user?.email}
              </span>
              <button
                onClick={logout}
                className={`flex items-center space-x-2 px-4 py-2 rounded-lg transition-colors ${
                  theme === 'dark'
                    ? 'bg-slate-700 text-slate-300 hover:bg-slate-600'
                    : 'bg-slate-100 text-slate-700 hover:bg-slate-200'
                }`}
              >
                <LogOut className="h-4 w-4" />
                <span>Logout</span>
              </button>
            </div>
          </div>
        </div>
      </header>

      <main className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        {/* Welcome Section */}
        <div className="mb-8">
          <h2 className={`text-3xl font-bold mb-2 ${theme === 'dark' ? 'text-white' : 'text-slate-900'}`}>
            Welcome back!
          </h2>
          <p className={`${theme === 'dark' ? 'text-slate-400' : 'text-slate-600'}`}>
            Select a module to begin your computational analysis
          </p>
        </div>

        {/* Module Grid */}
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 xl:grid-cols-4 gap-6">
          {modules.map((module) => {
            const Icon = module.icon;
            return (
              <Link key={module.name} to={module.path}>
                <div className={`module-card card hover:shadow-xl transition-all duration-300 cursor-pointer group ${
                  theme === 'dark' ? 'bg-slate-800 border-slate-700 hover:border-slate-600' : ''
                }`}>
                  <div className={`bg-gradient-to-br ${module.gradient} w-14 h-14 rounded-xl flex items-center justify-center mb-4 shadow-lg group-hover:scale-110 transition-transform duration-300`}>
                    <Icon className="h-7 w-7 text-white" />
                  </div>
                  <h3 className={`text-xl font-semibold mb-2 ${theme === 'dark' ? 'text-white' : 'text-slate-900'}`}>
                    {module.name}
                  </h3>
                  <p className={`text-sm ${theme === 'dark' ? 'text-slate-400' : 'text-slate-600'}`}>
                    {module.description}
                  </p>
                </div>
              </Link>
            );
          })}
        </div>

        {/* Stats Section */}
        <div className="mt-12 grid grid-cols-1 md:grid-cols-3 gap-6">
          <div className={`metric-card ${theme === 'dark' ? 'bg-slate-800 border border-slate-700' : 'bg-white border'} rounded-xl`}>
            <div className="absolute top-0 right-0 w-24 h-24 bg-blue-500 rounded-full -translate-y-1/2 translate-x-1/2 opacity-10" />
            <p className={`text-sm ${theme === 'dark' ? 'text-slate-400' : 'text-slate-600'} mb-1`}>Compute Quota</p>
            <p className={`text-3xl font-bold ${theme === 'dark' ? 'text-white' : 'text-slate-900'}`}>100 hours</p>
          </div>
          <div className={`metric-card ${theme === 'dark' ? 'bg-slate-800 border border-slate-700' : 'bg-white border'} rounded-xl`}>
            <div className="absolute top-0 right-0 w-24 h-24 bg-green-500 rounded-full -translate-y-1/2 translate-x-1/2 opacity-10" />
            <p className={`text-sm ${theme === 'dark' ? 'text-slate-400' : 'text-slate-600'} mb-1`}>Storage Quota</p>
            <p className={`text-3xl font-bold ${theme === 'dark' ? 'text-white' : 'text-slate-900'}`}>50 GB</p>
          </div>
          <div className={`metric-card ${theme === 'dark' ? 'bg-slate-800 border border-slate-700' : 'bg-white border'} rounded-xl`}>
            <div className="absolute top-0 right-0 w-24 h-24 bg-purple-500 rounded-full -translate-y-1/2 translate-x-1/2 opacity-10" />
            <p className={`text-sm ${theme === 'dark' ? 'text-slate-400' : 'text-slate-600'} mb-1`}>Account Status</p>
            <p className="text-3xl font-bold text-green-500">Active</p>
          </div>
        </div>

        {/* Quick Links / Recent */}
        <div className="mt-12">
          <h3 className={`text-xl font-semibold mb-4 ${theme === 'dark' ? 'text-white' : 'text-slate-900'}`}>
            Quick Access
          </h3>
          <div className={`${theme === 'dark' ? 'bg-slate-800 border-slate-700' : 'bg-white'} rounded-xl border p-6`}>
            <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
              <a 
                href="https://www.uniprot.org" 
                target="_blank" 
                rel="noopener noreferrer"
                className={`flex items-center space-x-3 p-3 rounded-lg transition-colors ${
                  theme === 'dark' ? 'hover:bg-slate-700' : 'hover:bg-slate-50'
                }`}
              >
                <div className="w-10 h-10 bg-blue-100 rounded-lg flex items-center justify-center">
                  <span className="text-blue-600 font-bold text-sm">U</span>
                </div>
                <div>
                  <p className={`font-medium ${theme === 'dark' ? 'text-white' : 'text-slate-900'}`}>UniProt</p>
                  <p className={`text-xs ${theme === 'dark' ? 'text-slate-400' : 'text-slate-500'}`}>Protein Database</p>
                </div>
              </a>
              <a 
                href="https://www.rcsb.org" 
                target="_blank" 
                rel="noopener noreferrer"
                className={`flex items-center space-x-3 p-3 rounded-lg transition-colors ${
                  theme === 'dark' ? 'hover:bg-slate-700' : 'hover:bg-slate-50'
                }`}
              >
                <div className="w-10 h-10 bg-green-100 rounded-lg flex items-center justify-center">
                  <span className="text-green-600 font-bold text-sm">PDB</span>
                </div>
                <div>
                  <p className={`font-medium ${theme === 'dark' ? 'text-white' : 'text-slate-900'}`}>RCSB PDB</p>
                  <p className={`text-xs ${theme === 'dark' ? 'text-slate-400' : 'text-slate-500'}`}>Structure Database</p>
                </div>
              </a>
              <a 
                href="https://alphafold.ebi.ac.uk" 
                target="_blank" 
                rel="noopener noreferrer"
                className={`flex items-center space-x-3 p-3 rounded-lg transition-colors ${
                  theme === 'dark' ? 'hover:bg-slate-700' : 'hover:bg-slate-50'
                }`}
              >
                <div className="w-10 h-10 bg-purple-100 rounded-lg flex items-center justify-center">
                  <span className="text-purple-600 font-bold text-sm">AF</span>
                </div>
                <div>
                  <p className={`font-medium ${theme === 'dark' ? 'text-white' : 'text-slate-900'}`}>AlphaFold</p>
                  <p className={`text-xs ${theme === 'dark' ? 'text-slate-400' : 'text-slate-500'}`}>AI Structures</p>
                </div>
              </a>
              <a 
                href="https://pubchem.ncbi.nlm.nih.gov" 
                target="_blank" 
                rel="noopener noreferrer"
                className={`flex items-center space-x-3 p-3 rounded-lg transition-colors ${
                  theme === 'dark' ? 'hover:bg-slate-700' : 'hover:bg-slate-50'
                }`}
              >
                <div className="w-10 h-10 bg-amber-100 rounded-lg flex items-center justify-center">
                  <span className="text-amber-600 font-bold text-sm">PC</span>
                </div>
                <div>
                  <p className={`font-medium ${theme === 'dark' ? 'text-white' : 'text-slate-900'}`}>PubChem</p>
                  <p className={`text-xs ${theme === 'dark' ? 'text-slate-400' : 'text-slate-500'}`}>Compound Database</p>
                </div>
              </a>
            </div>
          </div>
        </div>
      </main>

      {/* Footer */}
      <footer className={`mt-12 py-6 border-t ${theme === 'dark' ? 'bg-slate-800 border-slate-700' : 'bg-white'}`}>
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="flex items-center justify-between">
            <p className={`text-sm ${theme === 'dark' ? 'text-slate-400' : 'text-slate-500'}`}>
              © 2024 OpenBioBench. Built for the scientific community.
            </p>
            <div className="flex items-center space-x-4">
              <a href="#" className={`text-sm ${theme === 'dark' ? 'text-slate-400 hover:text-slate-300' : 'text-slate-500 hover:text-slate-700'}`}>
                Documentation
              </a>
              <a href="#" className={`text-sm ${theme === 'dark' ? 'text-slate-400 hover:text-slate-300' : 'text-slate-500 hover:text-slate-700'}`}>
                API Reference
              </a>
              <a href="https://github.com" target="_blank" rel="noopener noreferrer" className={`text-sm ${theme === 'dark' ? 'text-slate-400 hover:text-slate-300' : 'text-slate-500 hover:text-slate-700'}`}>
                GitHub
              </a>
            </div>
          </div>
        </div>
      </footer>
    </div>
  );
};
