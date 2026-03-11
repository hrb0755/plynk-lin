#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import resource
import subprocess
import time
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run one benchmark command and capture wall time plus peak RSS")
    parser.add_argument("--tool", required=True)
    parser.add_argument("--dataset-family", required=True)
    parser.add_argument("--dataset-name", required=True)
    parser.add_argument("--sample-count", required=True, type=int)
    parser.add_argument("--variant-count", required=True, type=int)
    parser.add_argument("--csv", required=True)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--shell-command", required=True)
    parser.add_argument("--poll-interval", type=float, default=0.05)
    return parser.parse_args()


def read_vm_rss_kb(pid: int) -> int | None:
    status_path = Path("/proc") / str(pid) / "status"
    try:
        for line in status_path.read_text(encoding="utf-8").splitlines():
            if line.startswith("VmRSS:"):
                parts = line.split()
                if len(parts) >= 2:
                    return int(parts[1])
    except FileNotFoundError:
        return None
    return None


def ensure_csv_header(path: Path) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "tool",
                "dataset_family",
                "dataset_name",
                "sample_count",
                "variant_count",
                "wall_seconds",
                "peak_rss_kb",
                "exit_code",
                "run_dir",
                "command",
            ]
        )


def append_csv_row(path: Path, row: list[object]) -> None:
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(row)


def main() -> int:
    ns = parse_args()
    run_dir = Path(ns.run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    csv_path = Path(ns.csv)
    ensure_csv_header(csv_path)

    stdout_path = run_dir / "stdout.log"
    stderr_path = run_dir / "stderr.log"
    shell_command = ns.shell_command

    before_rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    peak_polled = 0
    start = time.perf_counter()
    with stdout_path.open("w", encoding="utf-8", newline="") as stdout_handle, stderr_path.open(
        "w", encoding="utf-8", newline=""
    ) as stderr_handle:
        proc = subprocess.Popen(
            shell_command,
            shell=True,
            executable="/bin/bash",
            stdout=stdout_handle,
            stderr=stderr_handle,
        )
        while True:
            current = read_vm_rss_kb(proc.pid)
            if current is not None:
                peak_polled = max(peak_polled, current)
            if proc.poll() is not None:
                break
            time.sleep(ns.poll_interval)
        exit_code = proc.wait()
    wall_seconds = time.perf_counter() - start
    after_rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    peak_rss_kb = max(peak_polled, after_rss, before_rss) if os.name == "posix" else peak_polled

    row = [
        ns.tool,
        ns.dataset_family,
        ns.dataset_name,
        ns.sample_count,
        ns.variant_count,
        f"{wall_seconds:.6f}",
        peak_rss_kb,
        exit_code,
        str(run_dir),
        shell_command,
    ]
    append_csv_row(csv_path, row)

    print(
        f"{ns.tool}\t{ns.dataset_family}\t{ns.dataset_name}\t"
        f"samples={ns.sample_count}\tvariants={ns.variant_count}\t"
        f"wall_s={wall_seconds:.6f}\tpeak_rss_kb={peak_rss_kb}\texit_code={exit_code}"
    )
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
