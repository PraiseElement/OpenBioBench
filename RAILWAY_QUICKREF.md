# Railway Quick Reference Card

## 🚀 Quick Commands

```bash
# Install Railway CLI
npm install -g @railway/cli

# Login
railway login

# Initialize project
railway init

# Deploy
railway up

# View logs
railway logs

# Open dashboard
railway open
```

---

## 📋 Deployment Checklist

```
☐ 1. Install Railway CLI
☐ 2. Run deploy-railway.ps1
☐ 3. Set SECRET_KEY in backend
☐ 4. Set BACKEND_CORS_ORIGINS in backend
☐ 5. Set VITE_API_URL in frontend
☐ 6. Configure S3/R2 storage
☐ 7. Test backend health endpoint
☐ 8. Test frontend login
☐ 9. Update CORS with actual frontend URL
☐ 10. Set up custom domains (optional)
```

---

## 🔑 Essential Environment Variables

### Backend

```env
SECRET_KEY=<generate-32-char-hex>
BACKEND_CORS_ORIGINS=["https://your-frontend.railway.app"]
MINIO_ENDPOINT=s3.amazonaws.com
MINIO_ACCESS_KEY=<your-key>
MINIO_SECRET_KEY=<your-secret>
MINIO_BUCKET=openbiobench
```

### Frontend

```env
VITE_API_URL=https://your-backend.railway.app/api
```

---

## 🧪 Testing URLs

After deployment, test these:

- ✅ Frontend: `https://your-frontend.railway.app`
- ✅ Backend Health: `https://your-backend.railway.app/api/health`
- ✅ API Docs: `https://your-backend.railway.app/api/docs`

---

## 🛠️ Troubleshooting

| Issue                  | Solution                                     |
| ---------------------- | -------------------------------------------- |
| Backend won't start    | Check logs: `railway logs --service backend` |
| Frontend can't connect | Verify `VITE_API_URL` and CORS settings      |
| Database errors        | Check `DATABASE_URL` is set correctly        |
| File upload fails      | Verify S3/MinIO credentials                  |
| CORS errors            | Update `BACKEND_CORS_ORIGINS`                |

---

## 📞 Quick Links

- 📖 Full Guide: See `RAILWAY.md`
- 🔧 Env Template: See `.env.railway`
- 📄 Files List: See `RAILWAY_FILES.md`
- 🚂 Railway Docs: https://docs.railway.app
- 💬 Railway Discord: https://discord.gg/railway

---

## 💡 Pro Tips

1. Use `railway run` to execute commands in Railway environment
2. Set up GitHub auto-deploy: Railway → Settings → Deploy Triggers
3. Enable database backups in Railway dashboard
4. Use Railway's built-in metrics for monitoring
5. Start with 1 replica, scale up as needed

---

**Quick Start**: Run `.\deploy-railway.ps1` and follow the prompts! 🎉
