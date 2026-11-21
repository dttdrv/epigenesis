#!/usr/bin/env python3
"""Persistent stdio NDJSON compute kernel for the native Epigenesis desktop application.

Communicates strictly over private OS pipes (stdin/stdout) with zero open TCP sockets,
zero HTTP server, and zero browser dependencies. Exits cleanly on stdin EOF.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys
import traceback

REPO = Path(__file__).resolve().parents[1]
for p in (REPO, REPO / "studio"):
    sp = str(p)
    if sp not in sys.path:
        sys.path.insert(0, sp)

import epigenesis_studio  # noqa: E402
from epigenesis_studio import circuit, kreimer, pipeline_api, server, simulate  # noqa: E402
from epigenesis_studio.refusal import Refusal  # noqa: E402


def _handle(req: dict) -> dict:
    op = req.get("op")
    if op == "ping":
        return {"pong": True}
    if op == "info":
        info = server.model_info()
        info["repo"] = str(REPO)
        info["folder"] = str(server.FOLDER)
        info["files"] = server.folder_files()
        if pipeline_api.available():
            info["pipeline_presets"] = pipeline_api.presets_info()
        if simulate.available():
            info["site_table"] = simulate.site_table()
        if circuit.available():
            info["circuit_example"] = circuit.example_enhancer()
        info["example_edits"] = json.loads((server.HERE / "example-edits.json").read_text(encoding="utf-8"))
        return info
    if op == "set_folder":
        folder = Path(str(req["folder"])).expanduser().resolve()
        if not folder.is_dir():
            raise Refusal("no_file", f"folder not found: {folder}")
        server.FOLDER = folder
        return {"folder": str(server.FOLDER), "files": server.folder_files()}
    if op == "read_file":
        path = Path(str(req["path"])).expanduser().resolve()
        raw = path.read_bytes()
        if len(raw) > server.canonical.MAX_JSON_BYTES:
            raise Refusal("too_large", "file exceeds maximum supported size")
        text = raw.decode("utf-8")
        kind = server.classify(path.name, text) or "sequence"
        return {"name": path.name, "path": str(path), "text": text, "kind": kind}
    if op == "score":
        payload = req.get("request")
        if not isinstance(payload, dict):
            raise Refusal("bad_request", "missing request object")
        return server.score(payload)
    if op == "verify":
        receipt = req.get("receipt")
        if not isinstance(receipt, dict):
            raise Refusal("bad_request", "missing receipt object")
        return server.verify(receipt)
    raise Refusal("unsupported", f"unknown engine op: {op}")


def main() -> int:
    # Warm-load checkpoints once on startup
    try:
        server.predict.load()
        kreimer.load()
        if hasattr(server, "biophys") and server.biophys is not None:
            server.biophys.load()
    except Exception:
        pass

    # Signal readiness
    sys.stdout.write(json.dumps({"ready": True, "version": getattr(epigenesis_studio, "__version__", "0.3.0")}) + "\n")
    sys.stdout.flush()

    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        req_id = None
        try:
            req = json.loads(line)
            req_id = req.get("id")
            data = _handle(req)
            resp = {"id": req_id, "ok": True, "data": data}
        except Refusal as ref:
            resp = {"id": req_id, "ok": False, "code": ref.code, "error": str(ref), "params": ref.params}
        except Exception as err:
            resp = {"id": req_id, "ok": False, "code": "error", "error": str(err), "trace": traceback.format_exc(limit=4)}
        sys.stdout.write(json.dumps(resp, allow_nan=False) + "\n")
        sys.stdout.flush()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
