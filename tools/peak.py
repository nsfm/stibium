#!/usr/bin/env python3
"""Run a command; report wall time + peak RSS via /proc VmHWM."""
import subprocess, sys, time
t0 = time.time()
p = subprocess.Popen(sys.argv[1:])
peak = 0
while p.poll() is None:
    try:
        with open(f"/proc/{p.pid}/status") as f:
            for line in f:
                if line.startswith("VmHWM"):
                    peak = max(peak, int(line.split()[1]))
    except FileNotFoundError:
        break
    time.sleep(0.5)
wall = time.time() - t0
print(f"PEAK: wall {wall:.1f}s, max RSS {peak/1048576:.2f} GB")
sys.exit(p.returncode or 0)
