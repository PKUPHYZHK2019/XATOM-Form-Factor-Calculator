import os
import sys
import time
import sqlite3
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
    "Simulation/I_Dynamics/15fs_EventDyn/Cache/540eV/bound_transition_cache.sqlite"
)

DEFAULT_DECAY_CACHE_PATH = (
    "/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/"
    "Simulation/I_Dynamics/15fs_EventDyn/Cache/540eV/decay_lifetime_cache.sqlite"
)

DEFAULT_NONRES_CACHE_PATH = (
    "/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/"
    "Simulation/I_Dynamics/15fs_EventDyn/Cache/540eV/nonres_fpp_cache.sqlite"
)


# ============================================================
# Helper functions
# ============================================================

def load_unique_configs(events_sum_folder, photonE):
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
    unique_configs = [str(config) for config in unique_configs]

    return unique_configs


def load_finished_configs(output_file):
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
    if not os.path.exists(output_file):
        with open(output_file, "w") as f:
            f.write("config\tf0\tfp\tfpp\n")


def initialize_error_file(error_file):
    if not os.path.exists(error_file):
        with open(error_file, "w") as f:
            f.write("config\terror_type\terror_message\ttraceback\n")


def safe_one_line(text):
    return str(text).replace("\n", "\\n").replace("\t", "    ")


def calculate_formfac_with_retry(
    cfx,
    config,
    photonE,
    chunk_temp_storage_folder,
    bound_cache_path,
    decay_cache_path,
    nonres_cache_path,
    max_retries=10,
    wait_seconds=30,
):
    """
    Retry only temporary SQLite 'database is locked' errors.
    Other errors are raised immediately.
    """

    for attempt in range(max_retries + 1):
        try:
            return cfx.calculate_formfac(
                config,
                photonE,
                temp_storage_folder=chunk_temp_storage_folder,
                bound_cache_db_path=bound_cache_path,
                decay_cache_db_path=decay_cache_path,
                nonres_cache_db_path=nonres_cache_path,
            )

        except sqlite3.OperationalError as exc:
            message = str(exc).lower()

            if "database is locked" not in message:
                raise

            if attempt >= max_retries:
                raise

            sleep_time = wait_seconds * (attempt + 1)

            print(
                f"[SQL LOCK] database is locked. "
                f"Retry {attempt + 1}/{max_retries} after {sleep_time} s.",
                flush=True,
            )

            time.sleep(sleep_time)


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
        "--chunk_id",
        type=int,
        default=0,
        help="Chunk index."
    )

    parser.add_argument(
        "--n_chunks",
        type=int,
        default=10,
        help="Total number of chunks."
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
        help="Folder containing results_15fs_<photonE>eV.npz files."
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
    chunk_id = args.chunk_id
    n_chunks = args.n_chunks

    if n_chunks < 1:
        raise ValueError(f"n_chunks must be >= 1, got {n_chunks}")

    if chunk_id < 0 or chunk_id >= n_chunks:
        raise ValueError(
            f"chunk_id must satisfy 0 <= chunk_id < n_chunks. "
            f"Got chunk_id={chunk_id}, n_chunks={n_chunks}"
        )

    work_dir = args.work_dir
    config_files_folder = args.config_files_folder
    events_sum_folder = args.events_sum_folder

    bound_cache_path = args.bound_cache_path
    decay_cache_path = args.decay_cache_path
    nonres_cache_path = args.nonres_cache_path

    sys.path.append(work_dir)

    import calculate_formfac_xatom as cfx

    Path(bound_cache_path).parent.mkdir(parents=True, exist_ok=True)
    Path(decay_cache_path).parent.mkdir(parents=True, exist_ok=True)
    Path(nonres_cache_path).parent.mkdir(parents=True, exist_ok=True)

    # ============================================================
    # Permanent output folder
    # ============================================================

    photonE_folder = os.path.join(
        config_files_folder,
        f"{int(photonE)}eV"
    )

    Path(photonE_folder).mkdir(parents=True, exist_ok=True)

    output_file = os.path.join(
        photonE_folder,
        f"formfac_{int(photonE)}eV_chunk_{chunk_id:03d}.txt"
    )

    error_file = os.path.join(
        photonE_folder,
        f"formfac_{int(photonE)}eV_chunk_{chunk_id:03d}_errors.txt"
    )

    initialize_output_file(output_file)
    initialize_error_file(error_file)

    # ============================================================
    # Private temporary XATOM folder for this chunk
    # ============================================================

    chunk_temp_storage_folder = os.path.join(
        config_files_folder,
        f"{int(photonE)}eV",
        f"Temp_chunk_{int(photonE)}eV_{chunk_id:03d}"
    )

    Path(chunk_temp_storage_folder).mkdir(parents=True, exist_ok=True)

    # ============================================================
    # Load configs and split by chunk
    # ============================================================

    all_unique_configs = load_unique_configs(events_sum_folder, photonE)
    unique_configs = all_unique_configs[chunk_id::n_chunks]

    # Skip successful configs from this chunk output.
    finished_configs = load_finished_configs(output_file)

    # Also skip successful configs from old full non-chunk output.
    old_full_output_file = os.path.join(
        photonE_folder,
        f"formfac_{int(photonE)}eV.txt"
    )

    old_finished_configs = load_finished_configs(old_full_output_file)
    finished_configs = finished_configs.union(old_finished_configs)

    print("=" * 120, flush=True)
    print(f"Photon energy: {photonE} eV", flush=True)
    print(f"Chunk: {chunk_id}/{n_chunks}", flush=True)
    print(f"Total unique configs before chunking: {len(all_unique_configs)}", flush=True)
    print(f"Configs in this chunk: {len(unique_configs)}", flush=True)
    print(f"Already finished configs relevant to this chunk: {len(finished_configs)}", flush=True)
    print(f"Chunk temporary XATOM folder: {chunk_temp_storage_folder}", flush=True)
    print(f"Output file: {output_file}", flush=True)
    print(f"Error file: {error_file}", flush=True)
    print(f"Old full output file used for skipping: {old_full_output_file}", flush=True)
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

        print(
            f"[START] {i + 1}/{len(unique_configs)} | "
            f"photonE={photonE}, chunk={chunk_id}/{n_chunks}",
            flush=True,
        )

        try:
            f0, fp, fpp = calculate_formfac_with_retry(
                cfx=cfx,
                config=config,
                photonE=photonE,
                chunk_temp_storage_folder=chunk_temp_storage_folder,
                bound_cache_path=bound_cache_path,
                decay_cache_path=decay_cache_path,
                nonres_cache_path=nonres_cache_path,
            )

            with open(output_file, "a") as f:
                f.write(
                    f"{config}\t"
                    f"{float(f0):.12g}\t"
                    f"{float(fp):.12g}\t"
                    f"{float(fpp):.12g}\n"
                )

            n_success += 1
            finished_configs.add(config)

            print(
                f"[OK] {i + 1}/{len(unique_configs)} | "
                f"chunk={chunk_id}/{n_chunks}, "
                f"success={n_success}, failed={n_failed}, skipped={n_skipped}",
                flush=True,
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
                f"chunk={chunk_id}/{n_chunks}, "
                f"{type(exc).__name__}: {exc}",
                flush=True,
            )

    print("=" * 120, flush=True)
    print("Finished", flush=True)
    print(f"Photon energy: {photonE} eV", flush=True)
    print(f"Chunk: {chunk_id}/{n_chunks}", flush=True)
    print(f"Success: {n_success}", flush=True)
    print(f"Failed: {n_failed}", flush=True)
    print(f"Skipped existing: {n_skipped}", flush=True)
    print(f"Output file: {output_file}", flush=True)
    print(f"Error file: {error_file}", flush=True)
    print(f"Chunk temporary XATOM folder: {chunk_temp_storage_folder}", flush=True)
    print(f"Bound cache path: {bound_cache_path}", flush=True)
    print(f"Decay cache path: {decay_cache_path}", flush=True)
    print(f"Nonres cache path: {nonres_cache_path}", flush=True)
    print("=" * 120, flush=True)


if __name__ == "__main__":
    main()