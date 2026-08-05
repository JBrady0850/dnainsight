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

from flask import Flask, send_from_directory, jsonify
from backend import APP_VERSION
from backend.routes import api, MAX_UPLOAD_BYTES
from backend.database import init_db


def create_app() -> Flask:
    # Ensure the schema exists here rather than only in main(). init_db is
    # idempotent, and doing it at app-construction time means any WSGI host,
    # test client or embedded use gets a working database instead of failing on
    # the first query with "no such table: profiles".
    init_db()

    # v3.0 schema. Both are idempotent and additive: they CREATE TABLE IF NOT
    # EXISTS and never touch a table they did not create. Wrapped defensively
    # so that a v3 module missing from a partial checkout degrades to a working
    # v1 and v2 application rather than a server that will not boot.
    try:
        from backend.ledger import init_ledger
        init_ledger()
    except Exception as exc:  # pragma: no cover
        print(f"  WARNING: reclassification ledger unavailable ({exc}).")
    try:
        from backend.provenance import init_provenance
        init_provenance()
    except Exception as exc:  # pragma: no cover
        print(f"  WARNING: provenance store unavailable ({exc}).")

    app = Flask(__name__, static_folder=str(BASE_DIR / "frontend"))

    # Global request-size ceiling: Flask/Werkzeug aborts oversized requests
    # (413) before they are buffered, as a first line of defence.
    app.config["MAX_CONTENT_LENGTH"] = MAX_UPLOAD_BYTES

    @app.errorhandler(413)
    def _too_large(_e):
        mb = MAX_UPLOAD_BYTES // (1024 * 1024)
        return jsonify({"error": f"Upload exceeds the {mb} MB limit."}), 413

    # Register API blueprints. v1 first so its behaviour is authoritative for
    # every path it already owns; v2 adds the new surface described in
    # docs/API_V2.md without altering any existing endpoint.
    app.register_blueprint(api)
    try:
        from backend.routes_v2 import api_v2
        app.register_blueprint(api_v2)
    except Exception as exc:  # pragma: no cover
        # A missing v2 data file must degrade to a working v1 app, never a
        # server that will not boot.
        print(f"  WARNING: v2 API unavailable ({exc}). v1 endpoints still active.")
    try:
        from backend.routes_v3 import api_v3
        app.register_blueprint(api_v3)
    except Exception as exc:  # pragma: no cover
        # Same rule one level up. v3 adds ancestry, haplogroups, imputation and
        # the rest; none of it is required for a working v1 or v2 install.
        print(f"  WARNING: v3 API unavailable ({exc}). v1 and v2 still active.")

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

    # create_app() initialises the database, so no separate call is needed here.
    app = create_app()

    url = f"http://{args.host}:{args.port}"
    print(f"\n{'='*55}")
    print(f"  DNAInsight v{APP_VERSION}")
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
