"""
app.py -- DNAInsight application entry point.

Launches a local Flask web server and serves both the API and the
frontend single-page application. Automatically opens the browser.

Usage:
    python app.py          # default port 5050
    python app.py --port 8080
"""

import os
import sys
import argparse
import threading
import webbrowser
from pathlib import Path

# Ensure the project root is on sys.path regardless of CWD
BASE_DIR = Path(__file__).parent.resolve()
sys.path.insert(0, str(BASE_DIR))

from flask import Flask, send_from_directory
from backend.routes import api
from backend.database import init_db


def create_app() -> Flask:
    app = Flask(__name__, static_folder=str(BASE_DIR / "frontend"))

    # Register API blueprint
    app.register_blueprint(api)

    # Serve frontend SPA
    @app.route("/")
    def index():
        return send_from_directory(app.static_folder, "index.html")

    @app.route("/<path:path>")
    def static_files(path):
        return send_from_directory(app.static_folder, path)

    return app


def main():
    parser = argparse.ArgumentParser(description="DNAInsight — Personal DNA Analysis Tool")
    parser.add_argument("--port", type=int, default=5050, help="Port to listen on (default: 5050)")
    parser.add_argument("--no-browser", action="store_true", help="Do not open browser automatically")
    parser.add_argument("--host", default="127.0.0.1", help="Host to bind (default: 127.0.0.1)")
    args = parser.parse_args()

    # Initialize database
    init_db()

    app = create_app()

    url = f"http://{args.host}:{args.port}"
    print(f"\n{'='*55}")
    print(f"  DNAInsight v1.0")
    print(f"  Running at: {url}")
    print(f"  Press Ctrl+C to stop")
    print(f"{'='*55}\n")

    if not args.no_browser:
        def _open():
            import time
            time.sleep(1.2)
            webbrowser.open(url)
        threading.Thread(target=_open, daemon=True).start()

    app.run(host=args.host, port=args.port, debug=False, threaded=True)


if __name__ == "__main__":
    main()
