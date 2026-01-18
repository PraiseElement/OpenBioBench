/**
 * Theme Store - Manages application theme (light/dark)
 */
import { create } from 'zustand';
import { persist } from 'zustand/middleware';

interface ThemeState {
  theme: 'light' | 'dark';
  accentColor: string;
  fontSize: 'small' | 'medium' | 'large';
  useSystemTheme: boolean;
  setTheme: (theme: 'light' | 'dark') => void;
  toggleTheme: () => void;
  setAccentColor: (color: string) => void;
  setFontSize: (size: 'small' | 'medium' | 'large') => void;
  setUseSystemTheme: (use: boolean) => void;
}

export const useThemeStore = create<ThemeState>()(
  persist(
    (set, get) => ({
      theme: 'light',
      accentColor: '#3b82f6', // blue-500
      fontSize: 'medium',
      useSystemTheme: false,

      setTheme: (theme) => {
        set({ theme });
        updateDocumentTheme(theme);
      },

      toggleTheme: () => {
        const newTheme = get().theme === 'light' ? 'dark' : 'light';
        set({ theme: newTheme });
        updateDocumentTheme(newTheme);
      },

      setAccentColor: (color) => {
        set({ accentColor: color });
        document.documentElement.style.setProperty('--accent-color', color);
      },

      setFontSize: (fontSize) => {
        set({ fontSize });
        const sizes = { small: '14px', medium: '16px', large: '18px' };
        document.documentElement.style.setProperty('--base-font-size', sizes[fontSize]);
      },

      setUseSystemTheme: (use) => {
        set({ useSystemTheme: use });
        if (use) {
          const systemTheme = window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
          set({ theme: systemTheme });
          updateDocumentTheme(systemTheme);
        }
      },
    }),
    {
      name: 'openbiobench-theme',
      onRehydrateStorage: () => (state) => {
        if (state) {
          updateDocumentTheme(state.theme);
          if (state.accentColor) {
            document.documentElement.style.setProperty('--accent-color', state.accentColor);
          }
        }
      },
    }
  )
);

function updateDocumentTheme(theme: 'light' | 'dark') {
  if (theme === 'dark') {
    document.documentElement.classList.add('dark');
  } else {
    document.documentElement.classList.remove('dark');
  }
}

// Initialize system theme listener
if (typeof window !== 'undefined') {
  window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', (e) => {
    const store = useThemeStore.getState();
    if (store.useSystemTheme) {
      store.setTheme(e.matches ? 'dark' : 'light');
    }
  });
}
