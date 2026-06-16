import numpy as np
import sys
import os

photonE = 680
work_dir = '/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/Simulation/I_Dynamics/15fs_EventDyn'
sys.path.append('/path/to/the/example_file.py')
import calculate_formfac_xatom as cfx
import importlib
importlib.reload(cfx)
config_files_folder = '/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/Simulation/I_Dynamics/15fs_EventDyn/ConfigFiles'


events_sum_folder = f'/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/Simulation/I_Dynamics/15fs_EventDyn/Events_Sum'

# The result file has already considered the time limit of the pulse

f = os.path.join(events_sum_folder, f"results_15fs_{photonE}eV.npz")   # example

data = np.load(f, allow_pickle=True)

print(data.files)

for key in data.files:
    print("\n", key)
    print("shape =", np.shape(data[key]))
    print("dtype =", data[key].dtype)

    try:
        print("first few values:")
        print(data[key][:5])
    except:
        print(data[key])


unique_configs = np.unique(data['Current_Nodes'])
print(unique_configs.shape)



import os
import time
import json
import csv
import platform
import traceback
from pathlib import Path
from datetime import datetime


# ============================================================
# Benchmark configuration
# ============================================================

photonE = 690  # change this if needed

n_configs_to_test = 3
config_indices = [0, 1, 2]  # keep small; each run may take 20–60 s

benchmark_name = "calculate_formfac_baseline"

benchmark_output_folder = Path(config_files_folder) / "benchmarks"
benchmark_output_folder.mkdir(parents=True, exist_ok=True)

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

csv_path = benchmark_output_folder / f"{benchmark_name}_{timestamp}.csv"
json_path = benchmark_output_folder / f"{benchmark_name}_{timestamp}.json"
md_path = benchmark_output_folder / f"{benchmark_name}_{timestamp}.md"


# ============================================================
# Helper functions
# ============================================================

def short_config(config, max_len=80):
    """
    Shorten long electronic configuration strings for readable logs.
    """
    config = str(config)
    if len(config) <= max_len:
        return config
    return config[:max_len] + " ..."


def safe_float(x):
    """
    Convert values to float if possible.
    Useful if f0, f1, f2 are numpy scalars.
    """
    try:
        return float(x)
    except Exception:
        return str(x)


def run_single_benchmark(config, photonE, config_index, run_index):
    """
    Run one calculate_formfac benchmark and return a result dictionary.
    """

    temp_storage_folder = (
        Path(config_files_folder)
        / f"{photonE}eV"
        / "benchmark_runs"
        / f"config_{config_index}_run_{run_index}_{timestamp}"
    )
    temp_storage_folder.mkdir(parents=True, exist_ok=True)

    result = {
        "benchmark_name": benchmark_name,
        "timestamp": timestamp,
        "photonE_eV": photonE,
        "config_index": config_index,
        "run_index": run_index,
        "config_short": short_config(config),
        "temp_storage_folder": str(temp_storage_folder),
        "success": False,
        "runtime_s": None,
        "f0": None,
        "f1": None,
        "f2": None,
        "error_type": None,
        "error_message": None,
    }

    start = time.perf_counter()

    try:
        f0, f1, f2 = cfx.calculate_formfac(
            config,
            photonE,
            temp_storage_folder=str(temp_storage_folder),
        )

        end = time.perf_counter()

        result["success"] = True
        result["runtime_s"] = end - start
        result["f0"] = safe_float(f0)
        result["f1"] = safe_float(f1)
        result["f2"] = safe_float(f2)

    except Exception as e:
        end = time.perf_counter()

        result["success"] = False
        result["runtime_s"] = end - start
        result["error_type"] = type(e).__name__
        result["error_message"] = str(e)
        result["traceback"] = traceback.format_exc()

    return result


def summarize_results(results):
    """
    Compute summary statistics from successful runs.
    """

    successful = [r for r in results if r["success"]]
    failed = [r for r in results if not r["success"]]

    runtimes = [r["runtime_s"] for r in successful]

    summary = {
        "benchmark_name": benchmark_name,
        "timestamp": timestamp,
        "photonE_eV": photonE,
        "n_total": len(results),
        "n_success": len(successful),
        "n_failed": len(failed),
        "runtime_total_s": sum(r["runtime_s"] for r in results if r["runtime_s"] is not None),
        "runtime_mean_s": None,
        "runtime_min_s": None,
        "runtime_max_s": None,
    }

    if runtimes:
        summary["runtime_mean_s"] = sum(runtimes) / len(runtimes)
        summary["runtime_min_s"] = min(runtimes)
        summary["runtime_max_s"] = max(runtimes)

    return summary


def write_csv(results, csv_path):
    """
    Save benchmark results as CSV.
    """

    fieldnames = [
        "benchmark_name",
        "timestamp",
        "photonE_eV",
        "config_index",
        "run_index",
        "success",
        "runtime_s",
        "f0",
        "f1",
        "f2",
        "config_short",
        "temp_storage_folder",
        "error_type",
        "error_message",
    ]

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for r in results:
            row = {key: r.get(key, None) for key in fieldnames}
            writer.writerow(row)


def write_json(results, summary, json_path):
    """
    Save complete benchmark results as JSON.
    """

    output = {
        "summary": summary,
        "environment": {
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "processor": platform.processor(),
        },
        "results": results,
    }

    with open(json_path, "w") as f:
        json.dump(output, f, indent=2)


def make_markdown_report(results, summary):
    """
    Make a Markdown report that can be copied directly into a development log.
    """

    lines = []

    lines.append(f"## {datetime.now().strftime('%Y-%m-%d')} — [Benchmark] `{benchmark_name}`")
    lines.append("")
    lines.append("### Goal")
    lines.append("Establish baseline performance before modifying `calculate_formfac()`.")
    lines.append("")
    lines.append("### Benchmark setup")
    lines.append(f"- Photon energy: `{photonE} eV`")
    lines.append(f"- Number of configs tested: `{summary['n_total']}`")
    lines.append(f"- Successful runs: `{summary['n_success']}`")
    lines.append(f"- Failed runs: `{summary['n_failed']}`")
    lines.append(f"- Python version: `{platform.python_version()}`")
    lines.append(f"- Platform: `{platform.platform()}`")
    lines.append("")

    lines.append("### Runtime summary")
    if summary["runtime_mean_s"] is not None:
        lines.append(f"- Total runtime: `{summary['runtime_total_s']:.3f} s`")
        lines.append(f"- Mean runtime per successful config: `{summary['runtime_mean_s']:.3f} s`")
        lines.append(f"- Min runtime: `{summary['runtime_min_s']:.3f} s`")
        lines.append(f"- Max runtime: `{summary['runtime_max_s']:.3f} s`")
    else:
        lines.append("- No successful runs.")
    lines.append("")

    lines.append("### Per-config results")
    lines.append("")
    lines.append("| Config index | Success | Runtime / s | f0 | f1 | f2 |")
    lines.append("|---:|:---:|---:|---:|---:|---:|")

    for r in results:
        runtime_str = f"{r['runtime_s']:.3f}" if r["runtime_s"] is not None else "NA"
        f0_str = f"{r['f0']:.8g}" if isinstance(r["f0"], float) else str(r["f0"])
        f1_str = f"{r['f1']:.8g}" if isinstance(r["f1"], float) else str(r["f1"])
        f2_str = f"{r['f2']:.8g}" if isinstance(r["f2"], float) else str(r["f2"])

        lines.append(
            f"| {r['config_index']} | {r['success']} | {runtime_str} | "
            f"{f0_str} | {f1_str} | {f2_str} |"
        )

    lines.append("")
    lines.append("### Decision")
    lines.append("Baseline recorded. Future code changes should be compared against this benchmark.")
    lines.append("")
    lines.append("### Output files")
    lines.append(f"- CSV: `{csv_path}`")
    lines.append(f"- JSON: `{json_path}`")
    lines.append(f"- Markdown: `{md_path}`")

    return "\n".join(lines)


# ============================================================
# Run benchmark
# ============================================================

results = []

print("=" * 80)
print(f"Benchmark: {benchmark_name}")
print(f"Photon energy: {photonE} eV")
print(f"Config indices: {config_indices}")
print("=" * 80)

for run_index, config_index in enumerate(config_indices):
    config = unique_configs[config_index]

    print("")
    print(f"Running config_index = {config_index}")
    print(f"Config: {short_config(config)}")

    result = run_single_benchmark(
        config=config,
        photonE=photonE,
        config_index=config_index,
        run_index=run_index,
    )

    results.append(result)

    if result["success"]:
        print(
            f"Success | runtime = {result['runtime_s']:.3f} s | "
            f"f0 = {result['f0']:.8g}, "
            f"f1 = {result['f1']:.8g}, "
            f"f2 = {result['f2']:.8g}"
        )
    else:
        print(
            f"Failed | runtime = {result['runtime_s']:.3f} s | "
            f"{result['error_type']}: {result['error_message']}"
        )


# ============================================================
# Save results
# ============================================================

summary = summarize_results(results)

write_csv(results, csv_path)
write_json(results, summary, json_path)

markdown_report = make_markdown_report(results, summary)

with open(md_path, "w") as f:
    f.write(markdown_report)


# ============================================================
# Print final report
# ============================================================

print("")
print("=" * 80)
print("Benchmark finished")
print("=" * 80)
print("")
print(markdown_report)