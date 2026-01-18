# 🚀 Deploy OpenBioBench on Render

## Overview

Deploy your full-stack OpenBioBench application on Render with:

- ✅ **Free Tier** available
- ✅ **Automatic HTTPS**
- ✅ **Auto-deploys** from GitHub
- ✅ **PostgreSQL database** included
- ✅ **Global CDN**

---

## Quick Deploy (Using Blueprint)

### Step 1: Go to Render Dashboard

👉 [https://dashboard.render.com](https://dashboard.render.com)

### Step 2: Sign up/Login with GitHub

- Click **"Get Started"** → **"Login with GitHub"**
- Authorize Render to access your repositories

### Step 3: Deploy from Blueprint

1. Click **"New"** → **"Blueprint"**
2. Connect your **PraiseElement/OpenBioBench** repository
3. Render will automatically detect `render.yaml`
4. Click **"Apply"**
5. Wait for all services to deploy (5-10 minutes)

### Step 4: Access Your App

- **Frontend**: `https://openbiobench-frontend.onrender.com`
- **Backend API**: `https://openbiobench-backend.onrender.com/api/docs`

---

## Manual Setup (Alternative)

If you prefer manual setup:

### 1. Create PostgreSQL Database

1. Click **"New"** → **"PostgreSQL"**
2. Name: `openbiobench-db`
3. Database: `openbiobench`
4. User: `openbiobench_user`
5. Plan: **Free**
6. Click **"Create Database"**
7. Copy the **Internal Database URL**

### 2. Deploy Backend

1. Click **"New"** → **"Web Service"**
2. Connect to **PraiseElement/OpenBioBench**
3. Configure:
   - **Name**: `openbiobench-backend`
   - **Root Directory**: `backend`
   - **Runtime**: `Python 3`
   - **Build Command**: `pip install -r requirements.txt`
   - **Start Command**: `uvicorn app.main:app --host 0.0.0.0 --port $PORT`
   - **Plan**: **Free**

4. Add Environment Variables:

   ```
   DATABASE_URL = <paste Internal Database URL>
   SECRET_KEY = <generate-random-string>
   ACCESS_TOKEN_EXPIRE_MINUTES = 30
   CORS_ORIGINS = ["https://openbiobench-frontend.onrender.com"]
   ```

5. Click **"Create Web Service"**

### 3. Deploy Frontend

1. Click **"New"** → **"Web Service"**
2. Connect to **PraiseElement/OpenBioBench**
3. Configure:
   - **Name**: `openbiobench-frontend`
   - **Root Directory**: `frontend`
   - **Runtime**: `Node`
   - **Build Command**: `npm install && npm run build`
   - **Start Command**: `npm run preview -- --host --port $PORT`
   - **Plan**: **Free**

4. Add Environment Variables:

   ```
   VITE_API_URL = https://openbiobench-backend.onrender.com
   ```

5. Click **"Create Web Service"**

---

## Update Frontend API Configuration

Create `frontend/.env.production`:

```bash
VITE_API_URL=https://openbiobench-backend.onrender.com
```

Then update `frontend/src/services/api.ts`:

```typescript
const API_URL = import.meta.env.VITE_API_URL || "http://localhost:8000";
```

---

## Verify Deployment

1. **Check Backend Health**:

   ```
   https://openbiobench-backend.onrender.com/api/health
   ```

   Should return: `{"status": "healthy"}`

2. **Check Frontend**:
   Visit: `https://openbiobench-frontend.onrender.com`

3. **Test Registration**:
   - Create a new account
   - Login
   - Try one of the modules

---

## Important Notes

### Free Tier Limitations

- ⏸️ **Spins down after 15 minutes** of inactivity
- ⏰ **First request takes 30-60 seconds** to wake up
- 💾 **512 MB RAM** per service
- 🗄️ **1 GB PostgreSQL storage**

### Upgrade to Paid (Optional)

- 💰 **$7/month** per service
- ⚡ **Always on** (no spin down)
- 📈 **More RAM and storage**

---

## Troubleshooting

### Backend not starting?

Check logs: Dashboard → Backend Service → Logs

Common issues:

- Missing environment variables
- Database connection string incorrect
- Python version mismatch

### Frontend not connecting to backend?

- Check `CORS_ORIGINS` in backend includes frontend URL
- Verify `VITE_API_URL` is set correctly
- Check browser console for CORS errors

### Database connection failed?

- Use **Internal Database URL** (not External)
- Format: `postgresql://user:pass@host/dbname`

---

## Custom Domain (Optional)

1. Go to Frontend service
2. Click **"Settings"** → **"Custom Domain"**
3. Add your domain: `www.openbiobench.com`
4. Update DNS records as instructed
5. Update backend `CORS_ORIGINS`

---

## CI/CD (Auto-Deploy)

Render automatically deploys on every push to `master`:

```bash
git add .
git commit -m "Update feature"
git push origin master
```

Render will:

1. ✅ Build backend
2. ✅ Build frontend
3. ✅ Run migrations
4. ✅ Deploy both services

---

## Cost Estimate

| Service    | Plan       | Cost         |
| ---------- | ---------- | ------------ |
| Backend    | Free       | $0           |
| Frontend   | Free       | $0           |
| PostgreSQL | Free (1GB) | $0           |
| **Total**  |            | **$0/month** |

Or upgrade all to paid:

| Service    | Plan           | Cost          |
| ---------- | -------------- | ------------- |
| Backend    | Starter        | $7            |
| Frontend   | Starter        | $7            |
| PostgreSQL | Starter (10GB) | $7            |
| **Total**  |                | **$21/month** |

---

## Support

- 📖 [Render Docs](https://render.com/docs)
- 💬 [Render Community](https://community.render.com)
- 📧 Support: support@render.com

---

**Ready to deploy?** Just push your code and click "Apply" in Render! 🚀
