# 🧬 OpenBioBench

<div align="center">

![OpenBioBench Banner](https://img.shields.io/badge/OpenBioBench-Computational%20Biology%20Platform-blue?style=for-the-badge&logo=dna)

[![Python](https://img.shields.io/badge/Python-3.10+-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![React](https://img.shields.io/badge/React-18+-61DAFB?style=flat-square&logo=react&logoColor=black)](https://reactjs.org/)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.100+-009688?style=flat-square&logo=fastapi&logoColor=white)](https://fastapi.tiangolo.com/)
[![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)](LICENSE)

**A modern, full-stack web platform for computational biology and drug discovery research.**

[Features](#-features) • [Quick Start](#-quick-start) • [Documentation](#-documentation) • [Deployment](#-deployment)

</div>

---

## 🎯 Overview

OpenBioBench is an integrated computational biology workbench that provides researchers with essential tools for:

- 🔬 **Molecular Design** - Build and visualize drug-like molecules
- 📊 **ADMET Prediction** - Assess pharmacokinetic properties
- 🧬 **Sequence Analysis** - Multiple sequence alignment with conservation visualization
- 🎯 **Binding Site Detection** - Identify druggable pockets in proteins
- ⚗️ **Molecular Docking** - Predict protein-ligand binding modes
- 🌳 **Phylogenetic Analysis** - Construct evolutionary trees
- 🔧 **Structure Preparation** - Clean and prepare molecules for simulations

---

## ✨ Features

### 🌐 Database Integration

- **UniProt** - Search and fetch protein sequences
- **RCSB PDB** - Download crystal structures by ID
- **AlphaFold DB** - Access AI-predicted structures
- **PubChem** - Retrieve compound data and SMILES

### 🧪 Computational Modules

| Module                    | Description                                               |
| ------------------------- | --------------------------------------------------------- |
| **Ligand Builder**        | 2D/3D molecular design with RDKit                         |
| **ADMET Prediction**      | Absorption, Distribution, Metabolism, Excretion, Toxicity |
| **Sequence Alignment**    | Multiple sequence alignment with conservation scoring     |
| **Pocket Detection**      | Binding site identification and druggability analysis     |
| **Molecular Docking**     | Protein-ligand docking with pose prediction               |
| **Phylogenetics**         | Tree construction using NJ and UPGMA methods              |
| **Structure Preparation** | Protein/ligand cleaning and optimization                  |

### 🎨 User Experience

- Modern, responsive UI with Tailwind CSS
- Dark/Light theme support
- User authentication and project management
- File upload support (PDB, PDBQT, SDF, MOL2)

---

## 🚀 Quick Start

### Prerequisites

- **Python 3.10+**
- **Node.js 18+**
- **Docker & Docker Compose** (recommended)

### Option 1: Docker (Recommended)

```bash
# Clone the repository
git clone https://github.com/yourusername/OpenBioBench.git
cd OpenBioBench

# Start all services
docker-compose up -d

# Access the application
# Frontend: http://localhost:5173
# Backend API: http://localhost:8000/api/docs
```

### Option 2: Manual Setup

#### Backend Setup

```bash
cd backend

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Set up environment variables
cp .env.example .env
# Edit .env with your settings

# Run database migrations
alembic upgrade head

# Start the server
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```

#### Frontend Setup

```bash
cd frontend

# Install dependencies
npm install

# Start development server
npm run dev
```

#### Access the Application

- **Frontend**: http://localhost:5173
- **API Docs**: http://localhost:8000/api/docs
- **MinIO Console**: http://localhost:9001

---

## 📁 Project Structure

```
OpenBioBench/
├── backend/
│   ├── app/
│   │   ├── api/v1/          # API endpoints
│   │   ├── core/            # Config, security, database
│   │   ├── models/          # SQLAlchemy models
│   │   └── services/        # Business logic
│   ├── alembic/             # Database migrations
│   └── requirements.txt
├── frontend/
│   ├── src/
│   │   ├── modules/         # Feature modules
│   │   ├── components/      # Shared components
│   │   ├── store/           # State management
│   │   └── services/        # API clients
│   └── package.json
├── docker-compose.yml
└── README.md
```

---

## 🔧 Configuration

### Environment Variables

Create a `.env` file in the `backend/` directory:

```env
# Database
DATABASE_URL=postgresql+asyncpg://user:password@localhost:5432/openbiobench

# Security
SECRET_KEY=your-secret-key-here
ACCESS_TOKEN_EXPIRE_MINUTES=30

# CORS
CORS_ORIGINS=["http://localhost:5173","http://localhost:3000"]

# MinIO (Object Storage)
MINIO_ENDPOINT=localhost:9000
MINIO_ACCESS_KEY=minioadmin
MINIO_SECRET_KEY=minioadmin
```

---

## 🌍 Deployment

### No, You Don't Need Streamlit!

OpenBioBench is a **React + FastAPI** application. Streamlit is not required or recommended for this architecture.

### Deployment Options

| Platform             | Recommended For                        |
| -------------------- | -------------------------------------- |
| **Docker Compose**   | Self-hosted servers                    |
| **Railway**          | Easy cloud deployment                  |
| **Render**           | Free tier available                    |
| **AWS/GCP/Azure**    | Enterprise scale                       |
| **Vercel + Railway** | Frontend on Vercel, Backend on Railway |

### Railway Deployment (Recommended for Cloud)

OpenBioBench includes automated Railway configuration for seamless cloud deployment.

#### Quick Deploy

```bash
# Install Railway CLI
npm install -g @railway/cli

# Run deployment script
.\deploy-railway.ps1
```

#### Manual Setup

See [`RAILWAY.md`](RAILWAY.md) for detailed deployment instructions.

#### What Gets Deployed

- ✅ FastAPI backend with automatic migrations
- ✅ React frontend with production build
- ✅ PostgreSQL database (managed by Railway)
- ✅ Redis cache (managed by Railway)
- ✅ Automatic SSL certificates
- ✅ Health checks and auto-restart

### Docker Production Deployment

```bash
# Build frontend for production
cd frontend
npm run build

# Build Docker images
docker-compose -f docker-compose.prod.yml build

# Deploy
docker-compose -f docker-compose.prod.yml up -d
```

---

## 🧬 API Reference

### Authentication

```http
POST /api/v1/auth/register    # Create account
POST /api/v1/auth/login       # Get JWT token
```

### Database Access

```http
GET /api/v1/database/uniprot/{accession}   # Fetch UniProt sequence
GET /api/v1/database/pdb/{pdb_id}          # Fetch PDB structure
GET /api/v1/database/alphafold/{uniprot}   # Fetch AlphaFold structure
GET /api/v1/database/pubchem/{cid}         # Fetch compound data
```

### Computational Endpoints

```http
POST /api/v1/alignments/align        # Multiple sequence alignment
POST /api/v1/pockets/detect          # Binding pocket detection
POST /api/v1/docking/run             # Molecular docking
POST /api/v1/phylogeny/infer         # Phylogenetic tree
POST /api/v1/preparation/protein/clean   # Protein preparation
POST /api/v1/preparation/ligand/from_smiles   # Ligand preparation
```

Full API documentation available at `/api/docs` when the server is running.

---

## 🛠️ Tech Stack

### Backend

- **FastAPI** - Modern async Python web framework
- **SQLAlchemy** - Database ORM
- **PostgreSQL** - Primary database
- **RDKit** - Cheminformatics toolkit
- **BioPython** - Biological computation
- **MinIO** - Object storage

### Frontend

- **React 18** - UI framework
- **TypeScript** - Type safety
- **Tailwind CSS** - Styling
- **Zustand** - State management
- **Vite** - Build tool

---

## 🤝 Contributing

Contributions are welcome! Please read our [Contributing Guide](CONTRIBUTING.md) first.

```bash
# Fork the repository
# Create your feature branch
git checkout -b feature/amazing-feature

# Commit your changes
git commit -m 'Add amazing feature'

# Push to the branch
git push origin feature/amazing-feature

# Open a Pull Request
```

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- [RDKit](https://www.rdkit.org/) - Open-source cheminformatics
- [BioPython](https://biopython.org/) - Computational molecular biology
- [UniProt](https://www.uniprot.org/) - Protein sequence database
- [RCSB PDB](https://www.rcsb.org/) - Protein structure archive

---

<div align="center">

**Made with ❤️ for the scientific community**

[⬆ Back to Top](#-openbiobench)

</div>
