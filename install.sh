#!/usr/bin/env bash
# DNAInsight Installer -- macOS and Linux

set -e

echo ""
echo "============================================================"
echo "  DNAInsight Installer -- macOS / Linux"
echo "============================================================"
echo ""

# ============================================================
# STEP 0: Check / Auto-install Python 3
# ============================================================

if ! command -v python3 &>/dev/null; then
    echo "[WARN] Python 3 not found. Attempting automatic installation..."
    echo ""

    if command -v brew &>/dev/null; then
        echo "[AUTO] Homebrew detected. Installing Python 3..."
        brew install python3
    elif command -v apt &>/dev/null; then
        echo "[AUTO] apt detected. Installing Python 3..."
        sudo apt install -y python3 python3-pip python3-venv
    elif command -v apt-get &>/dev/null; then
        echo "[AUTO] apt-get detected. Installing Python 3..."
        sudo apt-get install -y python3 python3-pip python3-venv
    elif command -v dnf &>/dev/null; then
        echo "[AUTO] dnf detected. Installing Python 3..."
        sudo dnf install -y python3 python3-pip
    elif command -v yum &>/dev/null; then
        echo "[AUTO] yum detected. Installing Python 3..."
        sudo yum install -y python3 python3-pip
    elif command -v pacman &>/dev/null; then
        echo "[AUTO] pacman detected. Installing Python 3..."
        sudo pacman -S --noconfirm python python-pip
    else
        echo "[ERROR] Cannot auto-install Python 3. No supported package manager found."
        echo ""
        echo "  macOS:    brew install python3     (Homebrew: https://brew.sh)"
        echo "  Ubuntu:   sudo apt install python3 python3-pip"
        echo "  Fedora:   sudo dnf install python3 python3-pip"
        echo "  Arch:     sudo pacman -S python python-pip"
        echo ""
        echo "Then re-run this installer."
        exit 1
    fi

    if ! command -v python3 &>/dev/null; then
        echo "[ERROR] Python 3 install failed or not in PATH. Install manually and re-run."
        exit 1
    fi
fi

PYVER=$(python3 --version 2>&1)
echo "[OK] $PYVER detected."

# ============================================================
# STEP 1: Check / Auto-install pip
# ============================================================
echo ""
echo "[1/4] Checking pip..."

if ! python3 -m pip --version &>/dev/null; then
    echo "[AUTO] pip not found. Installing..."
    if command -v apt &>/dev/null || command -v apt-get &>/dev/null; then
        sudo apt install -y python3-pip 2>/dev/null || \
        sudo apt-get install -y python3-pip 2>/dev/null || \
        python3 -m ensurepip --upgrade
    elif command -v dnf &>/dev/null; then
        sudo dnf install -y python3-pip 2>/dev/null || python3 -m ensurepip --upgrade
    else
        python3 -m ensurepip --upgrade
    fi
fi

python3 -m pip install --upgrade pip --quiet 2>/dev/null || true
echo "[OK] pip ready."

# ============================================================
# STEP 2: Install Python dependencies
# ============================================================
echo ""
echo "[2/4] Installing Python dependencies..."

# Use --break-system-packages if needed (PEP 668 distros: Ubuntu 23.04+, Debian 12+)
if python3 -m pip install -r requirements.txt --quiet 2>/dev/null; then
    echo "[OK] Dependencies installed."
elif python3 -m pip install -r requirements.txt --quiet --break-system-packages 2>/dev/null; then
    echo "[OK] Dependencies installed (system pip)."
else
    echo "[AUTO] Creating virtual environment to avoid system pip conflict..."
    python3 -m venv .venv
    .venv/bin/pip install -r requirements.txt --quiet
    VENV_PYTHON=".venv/bin/python3"
    echo "[OK] Dependencies installed in virtual environment."
fi

# ============================================================
# STEP 3: Build bundled SNP reference
# ============================================================
echo ""
echo "[3/4] Building bundled SNP reference database..."

mkdir -p data uploads reports_output db

if [ -n "$VENV_PYTHON" ]; then
    $VENV_PYTHON data/build_reference.py
else
    python3 data/build_reference.py
fi
echo "[OK] SNP reference ready."

# ============================================================
# STEP 4: Create launcher script
# ============================================================
echo ""
echo "[4/4] Creating launch script..."

if [ -n "$VENV_PYTHON" ]; then
    cat > launch.sh << EOF
#!/usr/bin/env bash
cd "\$(dirname "\$0")"
.venv/bin/python3 app.py "\$@"
EOF
else
    cat > launch.sh << 'EOF'
#!/usr/bin/env bash
cd "$(dirname "$0")"
python3 app.py "$@"
EOF
fi

chmod +x launch.sh

echo ""
echo "============================================================"
echo "  Installation Complete!"
echo "============================================================"
echo ""
echo "To launch DNAInsight:"
echo "  ./launch.sh        (or)  python3 app.py"
echo ""
echo "Opens in your browser at http://127.0.0.1:5050"
echo ""
echo "TIP: Open DNAInsight and use Settings > Database to update"
echo "     ClinVar annotations monthly for the best results."
echo ""

read -r -p "Launch DNAInsight now? (y/N): " LAUNCH
if [[ "$LAUNCH" =~ ^[Yy]$ ]]; then
    if [ -n "$VENV_PYTHON" ]; then
        .venv/bin/python3 app.py &
    else
        python3 app.py &
    fi
    sleep 1.5
    if command -v open &>/dev/null; then
        open http://127.0.0.1:5050
    elif command -v xdg-open &>/dev/null; then
        xdg-open http://127.0.0.1:5050
    fi
    echo "Running. Press Ctrl+C to stop."
    wait
fi
