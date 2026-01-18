/**
 * Main App component with routing.
 */
import { BrowserRouter, Routes, Route, Navigate } from 'react-router-dom';
import { useEffect } from 'react';
import { useAuthStore } from './store/useAuthStore';
import { useThemeStore } from './store/useThemeStore';
import { Login } from './modules/auth/Login';
import { Register } from './modules/auth/Register';
import { Dashboard } from './modules/dashboard/Dashboard';
import { BuilderModule } from './modules/builder/BuilderModule';
import { ADMETModule } from './modules/admet/ADMETModule';
import { SequenceModule } from './modules/sequences/SequenceModule';
import { PocketModule } from './modules/pockets/PocketModule';
import { DockingModule } from './modules/docking/DockingModule';
import { PhylogenyModule } from './modules/phylogeny/PhylogenyModule';
import { PreparationModule } from './modules/preparation/PreparationModule';
import { SettingsModule } from './modules/settings/SettingsModule';
import './index.css';

function App() {
  const { loadUser, isAuthenticated } = useAuthStore();
  const { theme } = useThemeStore();

  useEffect(() => {
    loadUser();
  }, [loadUser]);

  // Apply theme class to document
  useEffect(() => {
    if (theme === 'dark') {
      document.documentElement.classList.add('dark');
    } else {
      document.documentElement.classList.remove('dark');
    }
  }, [theme]);

  return (
    <BrowserRouter>
      <Routes>
        <Route path="/login" element={<Login />} />
        <Route path="/register" element={<Register />} />
        
        <Route
          path="/"
          element={isAuthenticated ? <Dashboard /> : <Navigate to="/login" />}
        />
        <Route
          path="/builder"
          element={isAuthenticated ? <BuilderModule /> : <Navigate to="/login" />}
        />
        <Route
          path="/admet"
          element={isAuthenticated ? <ADMETModule /> : <Navigate to="/login" />}
        />
        <Route
          path="/sequences"
          element={isAuthenticated ? <SequenceModule /> : <Navigate to="/login" />}
        />
        <Route
          path="/pockets"
          element={isAuthenticated ? <PocketModule /> : <Navigate to="/login" />}
        />
        <Route
          path="/docking"
          element={isAuthenticated ? <DockingModule /> : <Navigate to="/login" />}
        />
        <Route
          path="/phylogeny"
          element={isAuthenticated ? <PhylogenyModule /> : <Navigate to="/login" />}
        />
        <Route
          path="/preparation"
          element={isAuthenticated ? <PreparationModule /> : <Navigate to="/login" />}
        />
        <Route
          path="/settings"
          element={isAuthenticated ? <SettingsModule /> : <Navigate to="/login" />}
        />
      </Routes>
    </BrowserRouter>
  );
}

export default App;



