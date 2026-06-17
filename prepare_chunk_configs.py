import os
from pathlib import Path

import numpy as np


WORK_DIR = (
    "/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/"
    "Simulation/I_Dynamics/15fs_EventDyn"
)

EVENTS_SUM_FOLDER = os.path.join(WORK_DIR, "Events_Sum")

CONFIG_CHUNKS_FOLDER = os.path.join(
    WORK_DIR,
    "Config_Chunks"
)

N_CHUNKS = 10

PHOTON_E_MIN = 540
PHOTON_E_MAX = 720


def main():
    Path(CONFIG_CHUNKS_FOLDER).mkdir(parents=True, exist_ok=True)

    for photonE in range(PHOTON_E_MAX, PHOTON_E_MIN - 1, -1):

        npz_file = os.path.join(
            EVENTS_SUM_FOLDER,
            f"results_15fs_{photonE}eV.npz"
        )

        if not os.path.exists(npz_file):
            print(f"{photonE} eV: missing {npz_file}, skip")
            continue

        out_dir = os.path.join(
            CONFIG_CHUNKS_FOLDER,
            f"{photonE}eV"
        )

        Path(out_dir).mkdir(parents=True, exist_ok=True)

        print("=" * 120)
        print(f"{photonE} eV: loading {npz_file}")

        data = np.load(npz_file, allow_pickle=True)

        if "Current_Nodes" not in data.files:
            raise KeyError(
                f"'Current_Nodes' not found in {npz_file}. "
                f"Available keys: {data.files}"
            )

        unique_configs = np.unique(data["Current_Nodes"])
        unique_configs = [str(config) for config in unique_configs]

        total_file = os.path.join(
            out_dir,
            f"total_configs_{photonE}eV.txt"
        )

        with open(total_file, "w") as f:
            f.write(f"{len(unique_configs)}\n")

        print(f"{photonE} eV: total unique configs = {len(unique_configs)}")

        for chunk_id in range(N_CHUNKS):
            chunk_configs = unique_configs[chunk_id::N_CHUNKS]

            chunk_file = os.path.join(
                out_dir,
                f"configs_{photonE}eV_chunk_{chunk_id:03d}.txt"
            )

            tmp_file = chunk_file + ".tmp"

            with open(tmp_file, "w") as f:
                for config in chunk_configs:
                    f.write(config + "\n")

            os.replace(tmp_file, chunk_file)

            print(
                f"{photonE} eV chunk {chunk_id:03d}: "
                f"{len(chunk_configs)} configs"
            )

    print("=" * 120)
    print("Done preparing config chunks")
    print("=" * 120)


if __name__ == "__main__":
    main()