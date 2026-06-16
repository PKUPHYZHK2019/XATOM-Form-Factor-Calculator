import os
import sys
import argparse
import traceback
from pathlib import Path

import numpy as np

# ============================================================
# User settings / default paths
# ============================================================

DEFAULT_WORK_DIR = (
    "/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/"
    "Simulation/I_Dynamics/15fs_EventDyn"
)

DEFAULT_CONFIG_FILES_FOLDER = (
    "/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/"
    "Simulation/I_Dynamics/15fs_EventDyn/ConfigFiles"
)

DEFAULT_EVENTS_SUM_FOLDER = (
    "/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/"
    "Simulation/I_Dynamics/15fs_EventDyn/Events_Sum"
)

DEFAULT_BOUND_CACHE_PATH = (
    "/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/"
    "Simulation/I_Dynamics/15fs_EventDyn/Cache/bound_transition_cache.sqlite"
)

DEFAULT_DECAY_CACHE_PATH = (
    "/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/"
    "Simulation/I_Dynamics/15fs_EventDyn/Cache/decay_lifetime_cache.sqlite"
)

DEFAULT_NONRES_CACHE_PATH = (
    "/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/"
    "Simulation/I_Dynamics/15fs_EventDyn/Cache/nonres_fpp_cache.sqlite"
)


# ============================================================
# Helper functions
# ============================================================

def load_unique_configs(events_sum_folder, photonE):
    """
    Load all unique electronic configurations for one photon energy.

    Expected input file:
        results_15fs_{photonE}eV.npz

    Expected key inside the npz:
        Current_Nodes
    """

    npz_file = os.path.join(
        events_sum_folder,
        f"results_15fs_{int(photonE)}eV.npz"
    )

    if not os.path.exists(npz_file):
        raise FileNotFoundError(f"Cannot find input file: {npz_file}")

    data = np.load(npz_file, allow_pickle=True)

    if "Current_Nodes" not in data.files:
        raise KeyError(
            f"'Current_Nodes' not found in {npz_file}. "
            f"Available keys: {data.files}"
        )

    unique_configs = np.unique(data["Current_Nodes"])

    # Convert possible numpy string/object types to normal Python strings
    unique_configs = [str(config) for config in unique_configs]

    return unique_configs


def load_finished_configs(output_file):
    """
    Read already-finished configs from the output file.

    This allows the job to resume if it was interrupted.
    """

    finished = set()

    if not os.path.exists(output_file):
        return finished

    with open(output_file, "r") as f:
        first_line = True

        for line in f:
            line = line.rstrip("\n")

            if first_line:
                first_line = False
                continue

            if not line.strip():
                continue

            config = line.split("\t")[0]
            finished.add(config)

    return finished


def initialize_output_file(output_file):
    """
    Create output file and write header if the file does not exist.
    """

    if not os.path.exists(output_file):
        with open(output_file, "w") as f:
            f.write("config\tf0\tfp\tfpp\n")


def initialize_error_file(error_file):
    """
    Create error file and write header if the file does not exist.
    """

    if not os.path.exists(error_file):
        with open(error_file, "w") as f:
            f.write("config\terror_type\terror_message\ttraceback\n")


def safe_one_line(text):
    """
    Convert error messages / tracebacks into one-line text for TSV output.
    """

    return str(text).replace("\n", "\\n").replace("\t", "    ")


# ============================================================
# Main simulation
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--photonE",
        type=float,
        required=True,
        help="Photon energy in eV."
    )

    parser.add_argument(
        "--work_dir",
        type=str,
        default=DEFAULT_WORK_DIR,
        help="Working directory containing calculate_formfac_xatom.py."
    )

    parser.add_argument(
        "--config_files_folder",
        type=str,
        default=DEFAULT_CONFIG_FILES_FOLDER,
        help="Folder for XATOM temporary files and output folders."
    )

    parser.add_argument(
        "--events_sum_folder",
        type=str,
        default=DEFAULT_EVENTS_SUM_FOLDER,
        help="Folder containing results_15fs_{photonE}eV.npz files."
    )

    parser.add_argument(
        "--bound_cache_path",
        type=str,
        default=DEFAULT_BOUND_CACHE_PATH,
        help="Path to bound-transition SQLite cache database."
    )

    parser.add_argument(
        "--decay_cache_path",
        type=str,
        default=DEFAULT_DECAY_CACHE_PATH,
        help="Path to decay-lifetime SQLite cache database."
    )

    parser.add_argument(
        "--nonres_cache_path",
        type=str,
        default=DEFAULT_NONRES_CACHE_PATH,
        help="Path to nonresonant fpp SQLite cache database."
    )

    args = parser.parse_args()

    photonE = args.photonE
    work_dir = args.work_dir
    config_files_folder = args.config_files_folder
    events_sum_folder = args.events_sum_folder

    bound_cache_path = args.bound_cache_path
    decay_cache_path = args.decay_cache_path
    nonres_cache_path = args.nonres_cache_path

    # Make sure Python can find calculate_formfac_xatom.py
    sys.path.append(work_dir)

    import calculate_formfac_xatom as cfx

    # Make sure SQL cache folders exist
    Path(bound_cache_path).parent.mkdir(parents=True, exist_ok=True)
    Path(decay_cache_path).parent.mkdir(parents=True, exist_ok=True)
    Path(nonres_cache_path).parent.mkdir(parents=True, exist_ok=True)

    # Output folder for this photon energy.
    # Permanent output files are saved in integer folders such as 540eV.
    photonE_folder = os.path.join(
        config_files_folder,
        f"{int(photonE)}eV"
    )

    Path(photonE_folder).mkdir(parents=True, exist_ok=True)

    output_file = os.path.join(
        photonE_folder,
        f"formfac_{int(photonE)}eV.txt"
    )

    error_file = os.path.join(
        photonE_folder,
        f"formfac_{int(photonE)}eV_errors.txt"
    )

    initialize_output_file(output_file)
    initialize_error_file(error_file)

    unique_configs = load_unique_configs(events_sum_folder, photonE)

    finished_configs = load_finished_configs(output_file)

    print("=" * 120, flush=True)
    print(f"Photon energy: {photonE} eV", flush=True)
    print(f"Number of unique configs: {len(unique_configs)}", flush=True)
    print(f"Already finished configs: {len(finished_configs)}", flush=True)
    print(f"Output file: {output_file}", flush=True)
    print(f"Error file: {error_file}", flush=True)
    print(f"Bound cache path: {bound_cache_path}", flush=True)
    print(f"Decay cache path: {decay_cache_path}", flush=True)
    print(f"Nonres cache path: {nonres_cache_path}", flush=True)
    print("=" * 120, flush=True)

    n_success = 0
    n_failed = 0
    n_skipped = 0

    for i, config in enumerate(unique_configs):
        if config in finished_configs:
            n_skipped += 1
            continue

        try:
            f0, fp, fpp = cfx.calculate_formfac(
                config,
                photonE,
                temp_storage_folder=config_files_folder,
                bound_cache_db_path=bound_cache_path,
                decay_cache_db_path=decay_cache_path,
                nonres_cache_db_path=nonres_cache_path,
            )

            with open(output_file, "a") as f:
                f.write(
                    f"{config}\t"
                    f"{float(f0):.12g}\t"
                    f"{float(fp):.12g}\t"
                    f"{float(fpp):.12g}\n"
                )

            n_success += 1

            print(
                f"[OK] {i + 1}/{len(unique_configs)} | "
                f"success={n_success}, failed={n_failed}, skipped={n_skipped}",
                flush=True
            )

        except Exception as exc:
            tb = traceback.format_exc()

            with open(error_file, "a") as f:
                f.write(
                    f"{config}\t"
                    f"{type(exc).__name__}\t"
                    f"{safe_one_line(exc)}\t"
                    f"{safe_one_line(tb)}\n"
                )

            n_failed += 1

            print(
                f"[FAILED] {i + 1}/{len(unique_configs)} | "
                f"{type(exc).__name__}: {exc}",
                flush=True
            )

    print("=" * 120, flush=True)
    print("Finished", flush=True)
    print(f"Photon energy: {photonE} eV", flush=True)
    print(f"Success: {n_success}", flush=True)
    print(f"Failed: {n_failed}", flush=True)
    print(f"Skipped existing: {n_skipped}", flush=True)
    print(f"Output file: {output_file}", flush=True)
    print(f"Error file: {error_file}", flush=True)
    print(f"Bound cache path: {bound_cache_path}", flush=True)
    print(f"Decay cache path: {decay_cache_path}", flush=True)
    print(f"Nonres cache path: {nonres_cache_path}", flush=True)
    print("=" * 120, flush=True)


if __name__ == "__main__":
    main()