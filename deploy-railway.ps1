#!/usr/bin/env pwsh
# Railway Deployment Script for OpenBioBench
# This script automates the Railway deployment process

Write-Host "🚂 OpenBioBench Railway Deployment Script" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan
Write-Host ""

# Check if Railway CLI is installed
Write-Host "Checking for Railway CLI..." -ForegroundColor Yellow
if (!(Get-Command railway -ErrorAction SilentlyContinue)) {
    Write-Host "❌ Railway CLI is not installed!" -ForegroundColor Red
    Write-Host ""
    Write-Host "Please install it first:" -ForegroundColor Yellow
    Write-Host "  npm install -g @railway/cli" -ForegroundColor White
    Write-Host ""
    Write-Host "Or visit: https://docs.railway.app/develop/cli" -ForegroundColor White
    exit 1
}
Write-Host "✅ Railway CLI found" -ForegroundColor Green
Write-Host ""

# Login to Railway
Write-Host "Logging into Railway..." -ForegroundColor Yellow
railway login

# Initialize project
Write-Host ""
Write-Host "Initializing Railway project..." -ForegroundColor Yellow
$projectExists = railway status 2>&1 | Select-String "Project:"

if (!$projectExists) {
    Write-Host "Creating new Railway project..." -ForegroundColor Yellow
    railway init
}
Write-Host "✅ Project initialized" -ForegroundColor Green
Write-Host ""

# Add PostgreSQL
Write-Host "Adding PostgreSQL database..." -ForegroundColor Yellow
railway add --plugin postgresql
Write-Host "✅ PostgreSQL added" -ForegroundColor Green
Write-Host ""

# Add Redis
Write-Host "Adding Redis cache..." -ForegroundColor Yellow
railway add --plugin redis
Write-Host "✅ Redis added" -ForegroundColor Green
Write-Host ""

# Setup Backend Service
Write-Host "Setting up Backend service..." -ForegroundColor Yellow
Set-Location backend
railway up --service backend
Set-Location ..
Write-Host "✅ Backend service deployed" -ForegroundColor Green
Write-Host ""

# Setup Frontend Service
Write-Host "Setting up Frontend service..." -ForegroundColor Yellow
Set-Location frontend
railway up --service frontend
Set-Location ..
Write-Host "✅ Frontend service deployed" -ForegroundColor Green
Write-Host ""

# Generate SECRET_KEY
Write-Host "Generating SECRET_KEY..." -ForegroundColor Yellow
$secretKey = -join ((48..57) + (65..90) + (97..122) | Get-Random -Count 32 | ForEach-Object {[char]$_})
Write-Host "Generated SECRET_KEY: $secretKey" -ForegroundColor Green
Write-Host ""

# Instructions for next steps
Write-Host "🎉 Deployment Complete!" -ForegroundColor Green
Write-Host ""
Write-Host "⚠️  IMPORTANT: You need to set these environment variables manually" -ForegroundColor Yellow
Write-Host "===============================================================" -ForegroundColor Yellow
Write-Host ""
Write-Host "Backend Service:" -ForegroundColor Cyan
Write-Host "  1. Go to Railway Dashboard → Backend Service → Variables" -ForegroundColor White
Write-Host "  2. Add the following variables:" -ForegroundColor White
Write-Host ""
Write-Host "     SECRET_KEY=$secretKey" -ForegroundColor White
Write-Host "     ALGORITHM=HS256" -ForegroundColor White
Write-Host "     ACCESS_TOKEN_EXPIRE_MINUTES=30" -ForegroundColor White
Write-Host "     DEBUG=False" -ForegroundColor White
Write-Host "     BACKEND_CORS_ORIGINS=[`"https://your-frontend.railway.app`"]" -ForegroundColor White
Write-Host ""
Write-Host "Frontend Service:" -ForegroundColor Cyan
Write-Host "  1. Go to Railway Dashboard → Frontend Service → Variables" -ForegroundColor White
Write-Host "  2. Add the following variable:" -ForegroundColor White
Write-Host ""
Write-Host "     VITE_API_URL=https://your-backend.railway.app/api" -ForegroundColor White
Write-Host ""
Write-Host "After setting variables:" -ForegroundColor Yellow
Write-Host "  1. Update BACKEND_CORS_ORIGINS with your actual frontend URL" -ForegroundColor White
Write-Host "  2. Redeploy both services" -ForegroundColor White
Write-Host ""
Write-Host "📚 For detailed instructions, see RAILWAY.md" -ForegroundColor Cyan
Write-Host ""
Write-Host "To view your deployment status:" -ForegroundColor Yellow
Write-Host "  railway status" -ForegroundColor White
Write-Host ""
Write-Host "To view logs:" -ForegroundColor Yellow
Write-Host "  railway logs --service backend" -ForegroundColor White
Write-Host "  railway logs --service frontend" -ForegroundColor White
Write-Host ""
