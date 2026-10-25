#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""run_all.py - executa script01/02/03 e grava as saidas em results/."""
import argparse
import contextlib
import hashlib
import io
import json
import platform
import runpy
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DEFAULT_SCRIPTS = ["script01.py", "script02.py", "script03.py"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def decode_result(value):
    """Aceita result = {"json": "<string JSON>"} ou um objeto JSON-serializavel."""
    if isinstance(value, dict) and "json" in value and len(value) == 1:
        inner = value["json"]
        if isinstance(inner, str):
            return json.loads(inner)
        return inner
    return value


def run_one(name, out_dir):
    path = ROOT / name
    info = {"script": name, "sha256": None, "status": "erro",
            "segundos": None, "arquivos": []}
    if not path.is_file():
        info["erro"] = "arquivo nao encontrado"
        return info
    info["sha256"] = sha256(path)
    stem = path.stem
    buf = io.StringIO()
    t0 = time.perf_counter()
    try:
        with contextlib.redirect_stdout(buf):
            # run_name != "__main__": os blocos `if __name__ == "__main__"`
            # (que imprimem o JSON) nao disparam; a saida vem de `result`.
            ns = runpy.run_path(str(path), run_name="run_all")
    except Exception:
        info["segundos"] = round(time.perf_counter() - t0, 3)
        info["erro"] = traceback.format_exc()
        return info
    info["segundos"] = round(time.perf_counter() - t0, 3)
    stdout_text = buf.getvalue()
    if "result" in ns:
        data = decode_result(ns["result"])
        target = out_dir / (stem + "_out.json")
        with open(target, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=1, default=str)
            f.write("\n")
        info["arquivos"].append(target.name)
    if stdout_text.strip():
        target = out_dir / (stem + "_stdout.txt")
        with open(target, "w", encoding="utf-8") as f:
            f.write(stdout_text)
        info["arquivos"].append(target.name)
    if not info["arquivos"]:
        info["erro"] = "o script nao definiu `result` nem imprimiu nada"
        return info
    info["status"] = "ok"
    return info


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Executa os scripts do article06 e grava results/.")
    parser.add_argument("scripts", nargs="*",
                        help="scripts a executar (padrao: script01..03)")
    parser.add_argument("--out", default="results",
                        help="pasta de saida (padrao: results)")
    args = parser.parse_args(argv)
    scripts = args.scripts or DEFAULT_SCRIPTS
    out_dir = Path(args.out)
    if not out_dir.is_absolute():
        out_dir = ROOT / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    report = []
    for name in scripts:
        print("==> %s ..." % name, flush=True)
        info = run_one(name, out_dir)
        report.append(info)
        if info["status"] == "ok":
            print("    ok (%.2f s): %s" % (info["segundos"], ", ".join(info["arquivos"])))
        else:
            msg = (info.get("erro") or "?").strip().splitlines()
            print("    FALHOU: %s" % (msg[-1] if msg else "?"))
    manifest = {"python": sys.version, "plataforma": platform.platform(),
                "scripts": report}
    with open(out_dir / "manifest.json", "w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=1)
        f.write("\n")
    return 0 if all(r["status"] == "ok" for r in report) else 1


if __name__ == "__main__":
    sys.exit(main())
