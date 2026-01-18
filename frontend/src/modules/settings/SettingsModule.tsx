/**
 * Settings Module.
 * User preferences, theme customization, and account settings.
 */
import React, { useState, useEffect } from 'react';
import { Link } from 'react-router-dom';
import { 
  Settings as SettingsIcon, 
  ArrowLeft, 
  Sun, 
  Moon, 
  User, 
  Palette, 
  Monitor,
  Save,
  Check,
  Bell,
  Shield,
  Database
} from 'lucide-react';
import { useThemeStore } from '../../store/useThemeStore';
import { useAuthStore } from '../../store/useAuthStore';

export const SettingsModule: React.FC = () => {
  const { theme, toggleTheme, accentColor, setAccentColor, fontSize, setFontSize, useSystemTheme, setUseSystemTheme } = useThemeStore();
  const { user } = useAuthStore();
  
  const [username, setUsername] = useState(user?.email?.split('@')[0] || '');
  const [email, setEmail] = useState(user?.email || '');
  const [notifications, setNotifications] = useState({
    jobComplete: true,
    weeklyDigest: false,
    updates: true,
  });
  const [saveSuccess, setSaveSuccess] = useState(false);

  const accentColors = [
    { name: 'Blue', value: '#3b82f6' },
    { name: 'Cyan', value: '#06b6d4' },
    { name: 'Emerald', value: '#10b981' },
    { name: 'Purple', value: '#8b5cf6' },
    { name: 'Rose', value: '#f43f5e' },
    { name: 'Orange', value: '#f97316' },
    { name: 'Amber', value: '#f59e0b' },
  ];

  const handleSave = () => {
    // In production, this would call the API
    setSaveSuccess(true);
    setTimeout(() => setSaveSuccess(false), 3000);
  };

  return (
    <div className={`min-h-screen ${theme === 'dark' ? 'bg-slate-900' : 'bg-slate-50'}`}>
      {/* Header */}
      <header className={`${theme === 'dark' ? 'bg-slate-800 border-slate-700' : 'bg-white'} border-b`}>
        <div className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 py-4">
          <div className="flex items-center space-x-4">
            <Link to="/" className={`${theme === 'dark' ? 'text-slate-400 hover:text-slate-200' : 'text-slate-600 hover:text-slate-900'}`}>
              <ArrowLeft className="h-6 w-6" />
            </Link>
            <div className="flex items-center space-x-3">
              <SettingsIcon className="h-8 w-8 text-blue-500" />
              <div>
                <h1 className={`text-2xl font-bold ${theme === 'dark' ? 'text-white' : 'text-slate-900'}`}>Settings</h1>
                <p className={`text-sm ${theme === 'dark' ? 'text-slate-400' : 'text-slate-600'}`}>Customize your experience</p>
              </div>
            </div>
          </div>
        </div>
      </header>

      <main className="max-w-4xl mx-auto px-4 sm:px-6 lg:px-8 py-8">
        <div className="space-y-8">
          {/* Profile Section */}
          <section className={`${theme === 'dark' ? 'bg-slate-800 border-slate-700' : 'bg-white'} rounded-xl border p-6`}>
            <div className="flex items-center space-x-3 mb-6">
              <User className="h-5 w-5 text-blue-500" />
              <h2 className={`text-lg font-semibold ${theme === 'dark' ? 'text-white' : 'text-slate-900'}`}>Profile</h2>
            </div>
            
            <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
              <div>
                <label className={`block text-sm font-medium mb-2 ${theme === 'dark' ? 'text-slate-300' : 'text-slate-700'}`}>
                  Display Name
                </label>
                <input
                  type="text"
                  value={username}
                  onChange={(e) => setUsername(e.target.value)}
                  className={`w-full px-4 py-2 rounded-lg border ${
                    theme === 'dark' 
                      ? 'bg-slate-700 border-slate-600 text-white' 
                      : 'bg-white border-slate-300'
                  } focus:ring-2 focus:ring-blue-500 focus:border-transparent`}
                />
              </div>
              <div>
                <label className={`block text-sm font-medium mb-2 ${theme === 'dark' ? 'text-slate-300' : 'text-slate-700'}`}>
                  Email Address
                </label>
                <input
                  type="email"
                  value={email}
                  onChange={(e) => setEmail(e.target.value)}
                  className={`w-full px-4 py-2 rounded-lg border ${
                    theme === 'dark' 
                      ? 'bg-slate-700 border-slate-600 text-white' 
                      : 'bg-white border-slate-300'
                  } focus:ring-2 focus:ring-blue-500 focus:border-transparent`}
                />
              </div>
            </div>
          </section>

          {/* Appearance Section */}
          <section className={`${theme === 'dark' ? 'bg-slate-800 border-slate-700' : 'bg-white'} rounded-xl border p-6`}>
            <div className="flex items-center space-x-3 mb-6">
              <Palette className="h-5 w-5 text-purple-500" />
              <h2 className={`text-lg font-semibold ${theme === 'dark' ? 'text-white' : 'text-slate-900'}`}>Appearance</h2>
            </div>

            {/* Theme Toggle */}
            <div className="mb-6">
              <label className={`block text-sm font-medium mb-3 ${theme === 'dark' ? 'text-slate-300' : 'text-slate-700'}`}>
                Theme
              </label>
              <div className="flex space-x-3">
                <button
                  onClick={() => { setUseSystemTheme(false); toggleTheme(); }}
                  className={`flex items-center space-x-2 px-4 py-2 rounded-lg border transition-all ${
                    theme === 'light' && !useSystemTheme
                      ? 'bg-blue-500 text-white border-blue-500'
                      : theme === 'dark'
                        ? 'bg-slate-700 border-slate-600 text-slate-300 hover:bg-slate-600'
                        : 'bg-slate-100 border-slate-300 hover:bg-slate-200'
                  }`}
                >
                  <Sun className="h-4 w-4" />
                  <span>Light</span>
                </button>
                <button
                  onClick={() => { setUseSystemTheme(false); if (theme === 'light') toggleTheme(); }}
                  className={`flex items-center space-x-2 px-4 py-2 rounded-lg border transition-all ${
                    theme === 'dark' && !useSystemTheme
                      ? 'bg-blue-500 text-white border-blue-500'
                      : theme === 'dark'
                        ? 'bg-slate-700 border-slate-600 text-slate-300 hover:bg-slate-600'
                        : 'bg-slate-100 border-slate-300 hover:bg-slate-200'
                  }`}
                >
                  <Moon className="h-4 w-4" />
                  <span>Dark</span>
                </button>
                <button
                  onClick={() => setUseSystemTheme(true)}
                  className={`flex items-center space-x-2 px-4 py-2 rounded-lg border transition-all ${
                    useSystemTheme
                      ? 'bg-blue-500 text-white border-blue-500'
                      : theme === 'dark'
                        ? 'bg-slate-700 border-slate-600 text-slate-300 hover:bg-slate-600'
                        : 'bg-slate-100 border-slate-300 hover:bg-slate-200'
                  }`}
                >
                  <Monitor className="h-4 w-4" />
                  <span>System</span>
                </button>
              </div>
            </div>

            {/* Accent Color */}
            <div className="mb-6">
              <label className={`block text-sm font-medium mb-3 ${theme === 'dark' ? 'text-slate-300' : 'text-slate-700'}`}>
                Accent Color
              </label>
              <div className="flex flex-wrap gap-3">
                {accentColors.map((color) => (
                  <button
                    key={color.value}
                    onClick={() => setAccentColor(color.value)}
                    className={`w-10 h-10 rounded-full flex items-center justify-center transition-transform hover:scale-110 ${
                      accentColor === color.value ? 'ring-2 ring-offset-2 ring-slate-400' : ''
                    }`}
                    style={{ backgroundColor: color.value }}
                    title={color.name}
                  >
                    {accentColor === color.value && <Check className="h-5 w-5 text-white" />}
                  </button>
                ))}
              </div>
            </div>

            {/* Font Size */}
            <div>
              <label className={`block text-sm font-medium mb-3 ${theme === 'dark' ? 'text-slate-300' : 'text-slate-700'}`}>
                Font Size
              </label>
              <div className="flex space-x-3">
                {(['small', 'medium', 'large'] as const).map((size) => (
                  <button
                    key={size}
                    onClick={() => setFontSize(size)}
                    className={`px-4 py-2 rounded-lg border capitalize transition-all ${
                      fontSize === size
                        ? 'bg-blue-500 text-white border-blue-500'
                        : theme === 'dark'
                          ? 'bg-slate-700 border-slate-600 text-slate-300 hover:bg-slate-600'
                          : 'bg-slate-100 border-slate-300 hover:bg-slate-200'
                    }`}
                  >
                    {size}
                  </button>
                ))}
              </div>
            </div>
          </section>

          {/* Notifications Section */}
          <section className={`${theme === 'dark' ? 'bg-slate-800 border-slate-700' : 'bg-white'} rounded-xl border p-6`}>
            <div className="flex items-center space-x-3 mb-6">
              <Bell className="h-5 w-5 text-amber-500" />
              <h2 className={`text-lg font-semibold ${theme === 'dark' ? 'text-white' : 'text-slate-900'}`}>Notifications</h2>
            </div>

            <div className="space-y-4">
              {[
                { key: 'jobComplete', label: 'Job completion alerts', description: 'Get notified when your computational jobs finish' },
                { key: 'weeklyDigest', label: 'Weekly digest', description: 'Receive a weekly summary of your activity' },
                { key: 'updates', label: 'Product updates', description: 'Stay informed about new features and improvements' },
              ].map((item) => (
                <div key={item.key} className="flex items-center justify-between py-3">
                  <div>
                    <p className={`font-medium ${theme === 'dark' ? 'text-white' : 'text-slate-900'}`}>{item.label}</p>
                    <p className={`text-sm ${theme === 'dark' ? 'text-slate-400' : 'text-slate-500'}`}>{item.description}</p>
                  </div>
                  <button
                    onClick={() => setNotifications(prev => ({ ...prev, [item.key]: !prev[item.key as keyof typeof prev] }))}
                    className={`relative w-12 h-6 rounded-full transition-colors ${
                      notifications[item.key as keyof typeof notifications] ? 'bg-blue-500' : theme === 'dark' ? 'bg-slate-600' : 'bg-slate-300'
                    }`}
                  >
                    <span
                      className={`absolute top-1 w-4 h-4 bg-white rounded-full transition-transform ${
                        notifications[item.key as keyof typeof notifications] ? 'translate-x-7' : 'translate-x-1'
                      }`}
                    />
                  </button>
                </div>
              ))}
            </div>
          </section>

          {/* Data & Privacy Section */}
          <section className={`${theme === 'dark' ? 'bg-slate-800 border-slate-700' : 'bg-white'} rounded-xl border p-6`}>
            <div className="flex items-center space-x-3 mb-6">
              <Shield className="h-5 w-5 text-green-500" />
              <h2 className={`text-lg font-semibold ${theme === 'dark' ? 'text-white' : 'text-slate-900'}`}>Data & Privacy</h2>
            </div>

            <div className="space-y-4">
              <div className={`p-4 rounded-lg ${theme === 'dark' ? 'bg-slate-700' : 'bg-slate-50'}`}>
                <div className="flex items-center justify-between">
                  <div className="flex items-center space-x-3">
                    <Database className="h-5 w-5 text-slate-400" />
                    <div>
                      <p className={`font-medium ${theme === 'dark' ? 'text-white' : 'text-slate-900'}`}>Storage Used</p>
                      <p className={`text-sm ${theme === 'dark' ? 'text-slate-400' : 'text-slate-500'}`}>2.4 GB of 50 GB</p>
                    </div>
                  </div>
                  <div className="w-32 h-2 bg-slate-200 rounded-full overflow-hidden">
                    <div className="w-[5%] h-full bg-blue-500 rounded-full" />
                  </div>
                </div>
              </div>

              <button className={`text-sm ${theme === 'dark' ? 'text-red-400 hover:text-red-300' : 'text-red-600 hover:text-red-700'}`}>
                Delete all my data
              </button>
            </div>
          </section>

          {/* Save Button */}
          <div className="flex justify-end">
            <button
              onClick={handleSave}
              className="flex items-center space-x-2 px-6 py-3 bg-blue-500 text-white rounded-lg hover:bg-blue-600 transition-colors"
            >
              {saveSuccess ? (
                <>
                  <Check className="h-5 w-5" />
                  <span>Saved!</span>
                </>
              ) : (
                <>
                  <Save className="h-5 w-5" />
                  <span>Save Changes</span>
                </>
              )}
            </button>
          </div>
        </div>
      </main>
    </div>
  );
};
