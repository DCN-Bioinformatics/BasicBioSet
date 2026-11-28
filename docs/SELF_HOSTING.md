# bioseq: Self-Hosting Guide

Hosting bioseq means serving a folder of static files. Any web server that can serve plain files will work, and you don't need Python, a database or a backend running on the server.

## Overview

bioseq is a small bioinformatics toolkit written in pure Python. It handles FASTA and FASTQ input, GC content, translation, ORFs, motifs, restriction sites, melting temperature and more. It has two front ends:

- **Terminal:** `main.py` prints reports and has an interactive menu.
- **Website:** `web/` runs the same Python code *inside the visitor's browser* using [Pyodide](https://pyodide.org). The analysis happens on the visitor's computer, so their sequences are never uploaded to your server.

The code is public at [github.com/DCN-Bioinformatics/BasicBioSet](https://github.com/DCN-Bioinformatics/BasicBioSet) and released under CC0 1.0 (public domain).

## Requirements

| What | Needed for | Notes |
| --- | --- | --- |
| A Linux server or VM | Hosting | Any distro. A small VM is enough because the server only sends files. |
| `git` | Downloading and updating | You can download a ZIP instead. |
| Python 3.10+ | Testing, terminal use, quick serving | Standard library only, nothing to `pip install`. |
| nginx or Apache | Production hosting | Recommended for a permanent public site. |
| Visitor internet access to `cdn.jsdelivr.net` | Website | The page loads the Pyodide engine (v0.26.4, about 10 MB on the first visit) from jsDelivr. |
| A modern browser | Website | Recent Chrome, Firefox, Edge or Safari. |

## Step 1: Download

Clone the repository into the folder you want to serve from. This guide uses `/srv/bioseq`.

```bash
sudo mkdir -p /srv/bioseq
sudo chown $USER /srv/bioseq
git clone https://github.com/DCN-Bioinformatics/BasicBioSet.git /srv/bioseq
cd /srv/bioseq
```

Without git, download the ZIP from the GitHub page (*Code → Download ZIP*) and extract it to the same place.

The parts that matter for hosting:

```
/srv/bioseq/
  web/        the website (index.html, app.js, style.css, bridge.py)
  bioseq/     the Python package the website loads into the browser
  data/       example sequences the website loads
  .nojekyll   only needed for GitHub Pages
```

**Serve the whole repository folder, not just `web/`.** The page fetches `../bioseq/*.py` and `../data/example.fasta`, so those folders have to be reachable next to `web/`.

## Step 2: Check that it works

From `/srv/bioseq`:

```bash
python3 -m unittest discover tests      # every test should pass
python3 main.py data/example.fasta      # prints a full report in the terminal
```

On a machine with a desktop, `python3 main.py --web` opens the website at `http://127.0.0.1:8000/web/`.

By default this server only listens on `127.0.0.1` (localhost), so **other computers can't reach it**. To open it to your network, add `--host 0.0.0.0` (see option A below). Step 3 covers serving to other people.

## Step 3: Serve it on your server

Use nginx (option B) for a permanent site. Option A is fine for a quick demo on a private network. Whichever option you use, the site lives at `/web/` under the folder you serve.

### Option A: Python's built-in server (quick test)

```bash
cd /srv/bioseq
python3 main.py --web --host 0.0.0.0 --port 8000
```

With a non-local `--host` it doesn't try to open a browser on the server. Without the bioseq code you can get the same with `python3 -m http.server 8000 --bind 0.0.0.0 --directory /srv/bioseq`.

Open `http://<server-ip>:8000/web/` from another computer. Press Ctrl+C to stop it. This server isn't hardened, so don't leave it exposed to the internet.

### Option B: nginx (recommended)

```bash
sudo apt install nginx            # Debian/Ubuntu; use dnf on Fedora/RHEL
```

Create `/etc/nginx/sites-available/bioseq`:

```nginx
server {
    listen 80;
    server_name bioseq.example.org;   # your domain, or _ for any

    root /srv/bioseq;

    location = / { return 302 /web/; }
    location ~ /\. { deny all; }      # hide .git, .github and other dotfiles

    location / {
        try_files $uri $uri/ =404;
    }
}
```

Enable it and reload:

```bash
sudo ln -s /etc/nginx/sites-available/bioseq /etc/nginx/sites-enabled/
sudo nginx -t && sudo systemctl reload nginx
```

The site is then at `http://bioseq.example.org/`, which redirects to `/web/`.

**Under a sub-path of an existing site:** all the page's paths are relative, so an alias works. Add this to your existing `server` block, and the site will be at `/bioseq/web/`:

```nginx
location /bioseq/ {
    alias /srv/bioseq/;
}
```

### Option C: Apache

```apache
<VirtualHost *:80>
    ServerName bioseq.example.org
    DocumentRoot /srv/bioseq
    RedirectMatch ^/$ /web/
    <Directory /srv/bioseq>
        Require all granted
    </Directory>
    <DirectoryMatch "/\.">
        Require all denied
    </DirectoryMatch>
</VirtualHost>
```

Save it as `/etc/apache2/sites-available/bioseq.conf`, then run `sudo a2ensite bioseq && sudo systemctl reload apache2`.

### Option D: GitHub Pages (no server at all)

Fork the repository, then go to *Settings → Pages → Deploy from a branch* and pick `main` and `/ (root)`. The site appears at `https://<user>.github.io/BasicBioSet/web/`. Keep the empty `.nojekyll` file, because without it GitHub skips `bioseq/__init__.py`.

## HTTPS, firewall and updates

**HTTPS:** for a public domain, get a free certificate from Let's Encrypt:

```bash
sudo apt install certbot python3-certbot-nginx   # or python3-certbot-apache
sudo certbot --nginx -d bioseq.example.org        # or --apache
```

**Firewall:** open the web ports, for example `sudo ufw allow 'Nginx Full'`. For option A, open port 8000 instead.

**Updating to a new version:**

```bash
cd /srv/bioseq && git pull
```

No restart is needed. Visitors' browsers may keep an old `web/app.js` for a few minutes, and that can briefly break the engine right after an update. To avoid this, add the block below inside the nginx `server` block so browsers always check for the newest files:

```nginx
location ~* \.(js|py|html)$ {
    add_header Cache-Control "no-cache";
}
```

## Troubleshooting

Most problems show up in the browser's developer console (F12 → Console and Network).

| Symptom | Likely cause | Fix |
| --- | --- | --- |
| Engine never starts and the page opens as `file:///...` | The page was opened from disk | Serve it over HTTP (Step 3). |
| 404 errors for `../bioseq/*.py` or `../data/example.fasta` | Only `web/` is being served | Point the server root at the whole repository folder. |
| 403 Forbidden | The web server can't read `/srv/bioseq` | `sudo chmod -R o+rX /srv/bioseq` |
| Can't connect from another computer | `main.py --web` was started without `--host 0.0.0.0`, or the firewall blocks the port | Add `--host 0.0.0.0` (or use option B or C) and open the port. |
| Engine fails to load on a closed network | Visitors can't reach `cdn.jsdelivr.net` | Allow that domain. Fully offline hosting isn't supported out of the box. |
| `ModuleNotFoundError` right after an update | The browser cached an old `app.js` | Hard-refresh (Ctrl+Shift+R) or add the no-cache rule above. |
| First load is slow | About 10 MB engine download | Normal. Later visits use the browser cache. |
| `bioseq/__init__.py` returns 404 on GitHub Pages | The `.nojekyll` file is missing | Restore the empty `.nojekyll` at the repository root. |
