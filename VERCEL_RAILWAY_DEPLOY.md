# 🚀 Deploy OpenBioBench: Vercel (Frontend) + Railway (Backend)

## ✅ Current Status

- **Frontend on Vercel**: ✅ Already deployed!
- **Backend on Railway**: 🔧 Need to configure

---

## Step 1: Deploy Backend on Railway

### A. Create Backend Service

1. Go to **Railway Dashboard**: https://railway.app/dashboard
2. Select your **OpenBioBench project**
3. Click **"+ New"** → **"GitHub Repo"**
4. Select **PraiseElement/OpenBioBench**

### B. Configure Backend Service

Click on the new service, then go to **Settings**:

**Root Directory:**

```
backend
```

**Config File Path:**

```
nixpacks.toml
```

**Start Command:** (if not auto-detected)

```
uvicorn app.main:app --host 0.0.0.0 --port $PORT
```

### C. Add Environment Variables

Go to **"Variables"** tab and add:

```
DATABASE_URL = ${{Postgres.DATABASE_URL}}
SECRET_KEY = your-secret-key-min-32-chars-long
ACCESS_TOKEN_EXPIRE_MINUTES = 30
CORS_ORIGINS = ["https://openbiobench.vercel.app"]
```

**Note:** Replace `openbiobench.vercel.app` with your actual Vercel URL.

To generate a secret key:

```bash
python -c "import secrets; print(secrets.token_urlsafe(32))"
```

---

## Step 2: Add PostgreSQL Database

1. In Railway project, click **"+ New"** → **"Database"** → **"PostgreSQL"**
2. Railway will auto-create `DATABASE_URL` variable
3. The backend service will automatically use it

---

## Step 3: Get Backend URL

After deployment completes:

1. Click on **backend service**
2. Go to **"Settings"**
3. Scroll to **"Networking"**
4. Click **"Generate Domain"**
5. Copy the URL (e.g., `openbiobench-backend-production.up.railway.app`)

---

## Step 4: Update Vercel Frontend

### A. Add Environment Variable in Vercel

1. Go to **Vercel Dashboard**: https://vercel.com
2. Select your **OpenBioBench** project
3. Go to **"Settings"** → **"Environment Variables"**
4. Add:
   - **Name**: `VITE_API_URL`
   - **Value**: `https://your-railway-backend-url.up.railway.app`
   - **Environments**: Production, Preview, Development

5. Click **"Save"**

### B. Redeploy Frontend

1. Go to **"Deployments"** tab
2. Click **"..."** on latest deployment
3. Click **"Redeploy"**

---

## Step 5: Update Backend CORS

Once you have the Vercel URL, update the backend CORS in Railway:

1. Go to Railway backend service
2. Go to **"Variables"** tab
3. Update `CORS_ORIGINS`:
   ```
   ["https://your-vercel-url.vercel.app"]
   ```
4. Redeploy backend

---

## ✅ Verification

Once both are deployed:

1. **Visit your Vercel URL**: `https://your-app.vercel.app`
2. **Register a new account**
3. **Login**
4. **Try a module** (e.g., Ligand Builder)

If everything works, you're done! 🎉

---

## 📊 Final Setup

| Service  | Platform           | URL                              | Cost                 |
| -------- | ------------------ | -------------------------------- | -------------------- |
| Frontend | Vercel             | `https://your-app.vercel.app`    | FREE                 |
| Backend  | Railway            | `https://backend.up.railway.app` | FREE ($5/mo credits) |
| Database | Railway PostgreSQL | Internal                         | FREE                 |

---

## 🐛 Troubleshooting

### CORS Error?

- Check `CORS_ORIGINS` in Railway backend includes your Vercel URL
- Redeploy backend after changing CORS

### API Connection Failed?

- Check `VITE_API_URL` in Vercel is correct
- Make sure Railway backend is deployed and running
- Verify backend health: `https://your-backend.up.railway.app/api/health`

### Backend Not Starting?

- Check Railway logs: Click service → "Deployments" → Latest → "View Logs"
- Verify `DATABASE_URL` is set
- Check `SECRET_KEY` is at least 32 characters

---

## 🎯 Next Steps

Your deployment URL will be:

- **App**: `https://openbiobench.vercel.app` (or your custom domain)
- **API Docs**: `https://your-backend.up.railway.app/api/docs`

Want to add a custom domain? I can help with that too!
