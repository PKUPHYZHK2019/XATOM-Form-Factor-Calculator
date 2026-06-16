# ============================================================
# Fundamental physical constants in SI units
# ============================================================

alpha = 7.2973525643e-3          # fine-structure constant, dimensionless

r_elec = 2.8179403262e-15            # classical electron radius, m

h = 6.62607015e-34               # Planck constant, J s

hbar = 1.054571817e-34           # reduced Planck constant, J s

c = 299792458.0                  # speed of light, m / s

e = 1.602176634e-19              # elementary charge, C

a0 = 5.29e-11                    # bohr radius



import subprocess
import shlex
from pathlib import Path
from collections import OrderedDict
import re
import os
import copy
import numpy as np
import pandas as pd
import sqlite3
import hashlib
import json
import io
import time



BOUND_CACHE_DB_PATH = (
    "/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/"
    "Simulation/I_Dynamics/15fs_EventDyn/Cache/bound_transition_cache.sqlite"
)

DECAY_CACHE_DB_PATH = (
    "/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/"
    "Simulation/I_Dynamics/15fs_EventDyn/Cache/decay_lifetime_cache.sqlite"
)

NONRES_CACHE_DB_PATH = (
    "/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/"
    "Simulation/I_Dynamics/15fs_EventDyn/Cache/nonres_fpp_cache.sqlite"
)
def has_nonzero_extended_orbitals(config, target_orbitals=("4f", "5d", "5f")):
    """
    Check whether selected orbitals have non-zero occupation in XATOM config.

    Example config segment:
        4f0,0
        5d0,0
        5f0,0

    Returns
    -------
    bool
        True if any target orbital has non-zero occupation.
    """

    parts = config.split("_")

    for part in parts:
        for orb in target_orbitals:
            if part.startswith(orb):
                occ_str = part[len(orb):]  # e.g. "0,0", "1,0", "0,2"

                occ_values = [
                    float(x)
                    for x in occ_str.split(",")
                    if x.strip() != ""
                ]

                if any(occ > 0 for occ in occ_values):
                    return True

    return False

def run_xatom(
    config,
    linewidth=1.0,
    PE_range="500-800",
    dE=1,
    element="I",
    output_dir="./ConfigFiles",
    is_resonance=True,
    verbose = False,
):
    """
    Run XATOM with external config, linewidth, and optional resonance mode.

    If 4f, 5d, or 5f has non-zero occupation, use:
        -rmax 100 -rmax_continuum 200

    Otherwise, do not add the rmax line.
    """

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Safe versions for shell command
    config_safe = shlex.quote(config)
    linewidth_safe = shlex.quote(str(linewidth))
    PE_range_safe = shlex.quote(str(PE_range))
    dE_safe = shlex.quote(str(dE))
    element_safe = shlex.quote(element)

    # Optional resonance line
    resonance_line = "    -resonance \\\n" if is_resonance else ""

    # Optional rmax line
    if has_nonzero_extended_orbitals(config):
        rmax_line = "    -rmax 100 -rmax_continuum 200 \\\n"
    else:
        rmax_line = ""

    # Optional label in filename
    resonance_label = "resonance" if is_resonance else "nonresonance"

    output_file = output_dir / f"{element}_{config}_{PE_range}_lw{linewidth}_{resonance_label}.out"
    output_file_safe = shlex.quote(str(output_file))



# source /opt/psi/Programming/anaconda/2019.07/conda/bin/activate
# conda activate /das/work/units/maloja/p19750/analysis/software/hankai-imageing

# module purge
# module load intel/22.2 
    cmd = f"""



/das/work/units/maloja/p21108/Hankai/Simulation/xraypac/xatom/src/xatom \\
    -s {element_safe} \\
    -config {config_safe} \\
{resonance_line}    -linewidth {linewidth_safe} \\
    -PE {PE_range_safe} \\
    -dE {dE_safe} \\
    -pcs \\
    -rel \\
{rmax_line}    -decay \\
    -o {output_file_safe}
"""

    if verbose:
        result = subprocess.run(
            ["bash", "-lc", cmd],
            capture_output=True,
            text=True
        )
    else:
        result = subprocess.run(
            ["bash", "-lc", cmd],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True
        )

    print(result.stdout)
    print(result.stderr)
    print("Return code:", result.returncode)
    print("Output file:", output_file)

    return output_file, result




# Atomic unit of time:
# 1 a.u. time = 2.4188843265857e-17 s = 0.024188843265857 fs
AU_TIME_FS = 0.024188843265857


def _read_text(file_path):
    """
    Read XATOM output file as text.
    """
    file_path = Path(file_path)
    return file_path.read_text(errors="ignore")


def _clean_xatom_line(line):
    """
    Remove leading comment marker and surrounding whitespace.
    This makes lines like:
        '#   P.E.(eV) Total ...'
    become:
        'P.E.(eV) Total ...'
    """
    return line.strip().lstrip("#").strip()
    
def _parse_lifetime_to_fs(lifetime_value, lifetime_unit):
    """
    Convert lifetime value with unit fs/as/ps/ns/s to fs.
    """
    lifetime_value = float(lifetime_value)
    unit = lifetime_unit.lower()

    if unit == "fs":
        return lifetime_value
    elif unit == "as":
        return lifetime_value * 1e-3
    elif unit == "ps":
        return lifetime_value * 1e3
    elif unit == "ns":
        return lifetime_value * 1e6
    elif unit == "s":
        return lifetime_value * 1e15
    else:
        raise ValueError(f"Unknown lifetime unit: {lifetime_unit}")


def extract_total_decay_info(file_path):
    """
    Extract total decay information from an XATOM output file.

    Returns
    -------
    dict
        Dictionary containing rates in a.u., lifetime in fs/as,
        and decay width in eV if available.
    """

    text = _read_text(file_path)

    info = {}

    # -------------------------------
    # Fluorescence rate
    # Example:
    # Total F rate =   8.85065E-06 a.u.  ( lifetime = 2733.0 fs )
    # -------------------------------
    m = re.search(
        r"Total\s+F\s+rate\s*=\s*([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s*a\.u\."
        r"(?:\s*\(\s*lifetime\s*=\s*([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s*([a-zA-Z]+)\s*\))?",
        text,
    )

    if m:
        rate = float(m.group(1))
        info["total_f_rate_au"] = rate

        if m.group(2) is not None:
            lifetime_fs = _parse_lifetime_to_fs(m.group(2), m.group(3))
        else:
            lifetime_fs = AU_TIME_FS / rate

        info["total_f_lifetime_fs"] = lifetime_fs
        info["total_f_lifetime_as"] = lifetime_fs * 1000

    # -------------------------------
    # Auger rate
    # Example:
    # Total A rate =   2.14478E-02 a.u.  ( lifetime = 1.1 fs )
    # -------------------------------
    m = re.search(
        r"Total\s+A\s+rate\s*=\s*([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s*a\.u\."
        r"(?:\s*\(\s*lifetime\s*=\s*([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s*([a-zA-Z]+)\s*\))?",
        text,
    )

    if m:
        rate = float(m.group(1))
        info["total_a_rate_au"] = rate

        if m.group(2) is not None:
            lifetime_fs = _parse_lifetime_to_fs(m.group(2), m.group(3))
        else:
            lifetime_fs = AU_TIME_FS / rate

        info["total_a_lifetime_fs"] = lifetime_fs
        info["total_a_lifetime_as"] = lifetime_fs * 1000

    # -------------------------------
    # Total decay rate
    # Example:
    # Total decay rate =   2.14567E-02 a.u.  ( lifetime = 1.1 fs )
    # -------------------------------
    m = re.search(
        r"Total\s+decay\s+rate\s*=\s*([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s*a\.u\."
        r"(?:\s*\(\s*lifetime\s*=\s*([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s*([a-zA-Z]+)\s*\))?",
        text,
    )

    if m:
        rate = float(m.group(1))
        info["total_decay_rate_au"] = rate

        # I recommend calculating lifetime from the rate,
        # because the printed lifetime may be rounded, e.g. "1.1 fs".
        lifetime_fs_from_rate = AU_TIME_FS / rate

        info["total_lifetime_fs"] = lifetime_fs_from_rate
        info["total_lifetime_as"] = lifetime_fs_from_rate * 1000

        if m.group(2) is not None:
            printed_lifetime_fs = _parse_lifetime_to_fs(m.group(2), m.group(3))
            info["printed_total_lifetime_fs"] = printed_lifetime_fs
            info["printed_total_lifetime_as"] = printed_lifetime_fs * 1000

    # -------------------------------
    # Total decay width
    # Example:
    # Total decay width =   5.83866E-01 eV
    # -------------------------------
    m = re.search(
        r"Total\s+decay\s+width\s*=\s*([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s*eV",
        text,
    )

    if m:
        info["total_decay_width_eV"] = float(m.group(1))

    return info




def extract_bound_bound_transitions(file_path):

    """

    Extract possible bound-to-bound transitions from XATOM output.

    Returns

    -------

    pandas.DataFrame

        Columns:

        - from_orb

        - to_orb

        - E_trans_eV

        - dipole_matrix_element

    """

    text = _read_text(file_path)

    start_marker = "Possible resonances for bound-to-bound transitions:"

    end_marker = "Photoabsorption cross section"

    start = text.find(start_marker)

    if start == -1:

        return pd.DataFrame(

            columns=["from_orb", "to_orb", "E_trans_eV", "dipole_matrix_element"]

        )

    end = text.find(end_marker, start)

    block = text[start:] if end == -1 else text[start:end]

    transitions = []

    # Allows optional leading "#"

    # Allows orbitals like 3d-, 3d+, 5s0, 4f, etc.

    pattern = re.compile(

        r"^\s*#?\s*"

        r"(?P<from_orb>\d+[a-zA-Z]+[0+\-]?)"

        r"\s+-\s+"

        r"(?P<to_orb>\d+[a-zA-Z]+[0+\-]?)"

        r"\s*:\s*"

        r"(?P<E_trans_eV>[+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)"

        r"\s+"

        r"(?P<dipole>[+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)"

    )

    for line in block.splitlines():

        m = pattern.match(line)

        if m:

            transitions.append(

                {

                    "from_orb": m.group("from_orb"),

                    "to_orb": m.group("to_orb"),

                    "E_trans_eV": float(m.group("E_trans_eV")),

                    "dipole_matrix_element": float(m.group("dipole")),

                }

            )

    return pd.DataFrame(

        transitions,

        columns=["from_orb", "to_orb", "E_trans_eV", "dipole_matrix_element"],

    )


def extract_pacs_at_energy(file_path, photon_energy, include_total=True):
    """
    Extract photoabsorption cross section contributions at a given photon energy.

    Parameters
    ----------
    file_path : str or pathlib.Path
        XATOM output file.
    photon_energy : int or float
        Photon energy in eV. Example: 540.
        It is assumed to exist in the file.
    include_total : bool
        If True, include "Total" in returned dict.
        If False, return only orbital contributions.

    Returns
    -------
    dict
        Example:
        {
            "Total": 0.834428,
            "4s0": 0.0783066,
            "3d-": 0.0,
            "3d+": 0.0,
            ...
        }
    """

    text = _read_text(file_path)

    start_marker = "Photoabsorption cross section (in Mb):"
    start = text.find(start_marker)

    if start == -1:
        raise ValueError("Could not find 'Photoabsorption cross section (in Mb):' block.")

    lines = text[start:].splitlines()

    labels = None
    target_energy = float(photon_energy)

    for line in lines:
        clean = _clean_xatom_line(line)

        if not clean:
            continue

        # Header line, possibly originally starting with "#"
        if clean.startswith("P.E.(eV)"):
            labels = clean.split()[1:]  # remove P.E.(eV)
            continue

        if labels is None:
            continue

        # Stop when table clearly ends
        if clean.startswith("Current date/time"):
            break

        parts = clean.split()

        # Data rows start with photon energy
        try:
            E = float(parts[0])
        except ValueError:
            continue

        if abs(E - target_energy) < 1e-6:
            values = [float(x) for x in parts[1:]]

            if len(values) != len(labels):
                raise ValueError(
                    f"Label/value mismatch at E={E}: "
                    f"{len(labels)} labels but {len(values)} values.\n"
                    f"labels = {labels}\n"
                    f"values = {values}"
                )

            result = dict(zip(labels, values))

            if not include_total:
                result.pop("Total", None)

            return result

    raise ValueError(f"Photon energy {photon_energy} eV not found in PACS block.")

import pandas as pd


def extract_pacs_table(file_path):
    """
    Extract photon-energy-dependent photoabsorption cross section table from XATOM output.

    Parameters
    ----------
    file_path : str or Path
        Path to XATOM .out file.

    Returns
    -------
    pacs_df : pandas.DataFrame
        Columns are:
        P.E.(eV), Total, 1s0, 2s0, 2p-, ...
        All values are floats.
    """

    with open(file_path, "r") as f:
        lines = f.readlines()

    # ------------------------------------------------------------
    # Find PACS block
    # ------------------------------------------------------------
    start_idx = None
    header_idx = None

    for i, line in enumerate(lines):
        if "Photoabsorption cross section" in line:
            start_idx = i
            break

    if start_idx is None:
        raise ValueError("Photoabsorption cross section block not found.")

    # The next line containing 'P.E.(eV)' is the column header
    for i in range(start_idx + 1, len(lines)):
        if "P.E.(eV)" in lines[i]:
            header_idx = i
            break

    if header_idx is None:
        raise ValueError("PACS header line with 'P.E.(eV)' not found.")

    # ------------------------------------------------------------
    # Parse header
    # Example:
    # #   P.E.(eV)        Total          1s0 ...
    # ------------------------------------------------------------
    header_line = lines[header_idx].strip()

    if header_line.startswith("#"):
        header_line = header_line[1:].strip()

    columns = header_line.split()

    # ------------------------------------------------------------
    # Parse numerical rows
    # Stop when reaching blank line, comment line, or non-numerical line
    # ------------------------------------------------------------
    rows = []

    for line in lines[header_idx + 1:]:
        stripped = line.strip()

        if stripped == "":
            break

        if stripped.startswith("#"):
            break

        parts = stripped.split()

        # Skip malformed lines
        if len(parts) != len(columns):
            break

        try:
            row = [float(x) for x in parts]
        except ValueError:
            break

        rows.append(row)

    if not rows:
        raise ValueError("PACS table was found, but no numerical rows were parsed.")

    pacs_df = pd.DataFrame(rows, columns=columns)

    return pacs_df
    
def _extract_rel_binding_energies_from_text(text, positive_binding=True, allow_empty=False):
    """
    Extract from block:
        Orbital energies with relativistic correction:
    """

    start_marker = "Orbital energies with relativistic correction:"
    start = text.find(start_marker)

    if start == -1:
        if allow_empty:
            return {}
        raise ValueError("Could not find relativistic orbital-energy block.")

    lines = text[start:].splitlines()

    binding_energies = {}
    in_table = False

    pattern = re.compile(
        r"^\s*"
        r"([0-9]+[spdfgh][0+-])"      # orbital, e.g. 3d-, 3d+, 5s0
        r"\s+"
        r"([0-9]+)"                   # n_occ
        r"\s+"
        r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)"  # E(eV)
    )

    for line in lines:
        clean = _clean_xatom_line(line)

        if not clean:
            if in_table and binding_energies:
                break
            continue

        if clean.startswith("n_occ"):
            in_table = True
            continue

        if not in_table:
            continue

        # Stop after the relativistic table
        if clean.startswith("SCF_REL") or clean.startswith("E(0)") or clean.startswith("E_TOT_REL"):
            break

        m = pattern.match(clean)

        if m:
            orb = m.group(1)
            E_eV = float(m.group(3))
            binding_energies[orb] = abs(E_eV) if positive_binding else E_eV

    if not binding_energies and not allow_empty:
        raise ValueError("Relativistic orbital-energy table was found, but no entries were parsed.")

    return binding_energies


def _extract_nonrel_binding_energies_from_text(text, positive_binding=True, allow_empty=False):
    """
    Extract from block:
        Orbital energies:
    """

    start_marker = "Orbital energies:"
    start = text.find(start_marker)

    if start == -1:
        if allow_empty:
            return {}
        raise ValueError("Could not find non-relativistic orbital-energy block.")

    lines = text[start:].splitlines()

    binding_energies = {}
    in_table = False

    pattern = re.compile(
        r"^\s*"
        r"([0-9]+[spdfgh])"           # orbital, e.g. 1s, 2p, 3d
        r"\s+"
        r"([0-9]+)"                   # n_occ
        r"\s+"
        r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)"  # E(a.u.)
        r"\s+"
        r"([+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)"  # E(eV)
    )

    for line in lines:
        clean = _clean_xatom_line(line)

        if not clean:
            if in_table and binding_energies:
                break
            continue

        if clean.startswith("n_occ"):
            in_table = True
            continue

        if not in_table:
            continue

        # Stop after the nonrelativistic table
        if clean.startswith("SCF") or clean.startswith("Energies with relativistic correction"):
            break

        m = pattern.match(clean)

        if m:
            orb = m.group(1)
            E_eV = float(m.group(4))
            binding_energies[orb] = abs(E_eV) if positive_binding else E_eV

    if not binding_energies and not allow_empty:
        raise ValueError("Non-relativistic orbital-energy table was found, but no entries were parsed.")

    return binding_energies

def extract_binding_energies(file_path, positive_binding=True, prefer_relativistic=True):
    """
    Extract orbital binding energies from XATOM output.

    Parameters
    ----------
    file_path : str or pathlib.Path
        XATOM output file.
    positive_binding : bool
        If True, return positive binding energies.
        Example:
            {"3d-": 661.61}
        If False, return raw orbital energies.
        Example:
            {"3d-": -661.61}
    prefer_relativistic : bool
        If True, first try to read the relativistic orbital-energy block.
        If not found or empty, fall back to the non-relativistic block.

    Returns
    -------
    dict
        Relativistic example:
        {
            "1s0": 33237.30,
            "2p-": 4898.72,
            "2p+": 4622.52,
            "3d-": 661.61,
            "3d+": 649.33,
            ...
        }

        Non-relativistic example:
        {
            "1s": 31840.21,
            "2p": 4556.40,
            "3d": 643.78,
            ...
        }
    """

    text = _read_text(file_path)

    if prefer_relativistic:
        rel = _extract_rel_binding_energies_from_text(
            text,
            positive_binding=positive_binding,
            allow_empty=True,
        )

        if len(rel) > 0:
            return rel

    nonrel = _extract_nonrel_binding_energies_from_text(
        text,
        positive_binding=positive_binding,
        allow_empty=True,
    )

    if len(nonrel) > 0:
        return nonrel

    raise ValueError("Could not extract either relativistic or non-relativistic binding energies.")


def config2dict(config_str):
    """
    Convert XATOM configuration string to an OrderedDict with split relativistic orbitals.

    Example
    -------
    "1s2_2s2_2p2,4_3d4,5"

    becomes

    OrderedDict({
        "1s0": 2,
        "2s0": 2,
        "2p-": 2,
        "2p+": 4,
        "3d-": 4,
        "3d+": 5,
    })

    Notes
    -----
    s orbitals are stored as "s0".
    p/d/f/... orbitals with comma occupations are stored as "-" and "+".
    """

    config_dict = OrderedDict()

    parts = config_str.strip().split("_")

    for part in parts:
        part = part.strip()
        if not part:
            continue

        m = re.match(r"^(\d+)([a-zA-Z]+)(.*)$", part)
        if m is None:
            raise ValueError(f"Cannot parse config part: {part}")

        n = m.group(1)
        l = m.group(2)
        occ_str = m.group(3)

        if "," in occ_str:
            occs = [int(x) for x in occ_str.split(",")]

            if len(occs) != 2:
                raise ValueError(f"Expected two split occupations in: {part}")

            config_dict[f"{n}{l}-"] = occs[0]
            config_dict[f"{n}{l}+"] = occs[1]

        else:
            config_dict[f"{n}{l}0"] = int(occ_str)

    return config_dict


def dict2config(config_dict):
    """
    Convert split relativistic orbital dict back to XATOM configuration string.

    Example
    -------
    OrderedDict({
        "1s0": 2,
        "2s0": 2,
        "2p-": 2,
        "2p+": 4,
        "3d-": 4,
        "3d+": 5,
    })

    becomes

    "1s2_2s2_2p2,4_3d4,5"
    """

    parts = []
    used = set()

    orbital_pattern = re.compile(r"^(\d+)([a-zA-Z]+)([0+\-])$")

    keys = list(config_dict.keys())

    for key in keys:
        if key in used:
            continue

        m = orbital_pattern.match(key)
        if m is None:
            raise ValueError(f"Invalid split orbital key: {key}")

        n, l, split = m.groups()

        if split == "0":
            parts.append(f"{n}{l}{config_dict[key]}")
            used.add(key)

        elif split == "-":
            plus_key = f"{n}{l}+"

            if plus_key not in config_dict:
                raise ValueError(f"Missing matching plus orbital for {key}: expected {plus_key}")

            parts.append(f"{n}{l}{config_dict[key]},{config_dict[plus_key]}")
            used.add(key)
            used.add(plus_key)

        elif split == "+":
            minus_key = f"{n}{l}-"

            if minus_key not in config_dict:
                raise ValueError(f"Missing matching minus orbital for {key}: expected {minus_key}")

            # If '+' appears before '-', still output the pair once.
            parts.append(f"{n}{l}{config_dict[minus_key]},{config_dict[key]}")
            used.add(minus_key)
            used.add(key)

    return "_".join(parts)




def kk_fp_single_from_fpp(PhotonEs, fpp_values, E_eval):

    """
    Calculate f'(E_eval) from f''(E') using the Kramers-Kronig relation:

        f'(E) = 2/pi * P ∫ E' f''(E') / (E^2 - E'^2) dE'

    Important sign convention
    -------------------------
    In this code, `fpp_values` is allowed to be the signed imaginary part
    used in the Kramers-Heisenberg convention.

    In the later ionization calculation, you define:

        fpp_ionization_total = -PACS * E / 69.9

    Therefore absorption corresponds to negative f''.

    The final line of this function returns:

        return -float(fp_value)

    This extra minus sign converts the standard positive-absorption
    Kramers-Kronig result into the sign convention used by the rest
    of this code.

    Principal value treatment
    -------------------------
    The integrand contains a formal singularity when E' = E:

        denominator = E^2 - E'^2 = 0

    Instead of directly dividing by zero, the code uses singularity
    subtraction:

        g(E') = E' f''(E')

        g(E') / (E^2 - E'^2)
        =
        [g(E') - g(E)] / (E^2 - E'^2)
        +
        g(E) / (E^2 - E'^2)

    The first term is finite at E' = E.
    The second term is integrated analytically as a principal-value
    integral over the finite energy range [E_min, E_max].

    Parameters
    ----------
    PhotonEs : array-like
        Photon energy grid in eV.
        This grid can be non-uniform, but it should cover the energy
        region where f'' contributes significantly.

    fpp_values : array-like
        f'' values on the same energy grid.
        In this code, these values are usually signed negative values
        for absorption.

    E_eval : float
        The photon energy where f' should be evaluated.

    Returns
    -------
    fp_value : float
        KK-calculated f'(E_eval), with the sign convention used in this code.
    """

    # Convert input energy grid and f'' values into NumPy arrays.
    # This ensures that vectorized operations such as multiplication,
    # sorting, masking, and interpolation work consistently.
    PhotonEs = np.asarray(PhotonEs, dtype=float)
    fpp_values = np.asarray(fpp_values, dtype=float)

    # Make sure the evaluation energy is a scalar float.
    E_eval = float(E_eval)

    # The KK integral requires one f'' value for each photon energy.
    # If their shapes do not match, the integral is ill-defined.
    if PhotonEs.shape != fpp_values.shape:
        raise ValueError("PhotonEs and fpp_values must have the same shape.")
        
    # Sort energy grid.
    # np.trapezoid(integrand, E_grid) assumes the x-axis is ordered.
    # Non-uniform spacing is okay, but unsorted spacing would give wrong signs
    # and wrong interval widths.
    sort_idx = np.argsort(PhotonEs)
    E_grid = PhotonEs[sort_idx]
    fpp_grid = fpp_values[sort_idx]
    
    # Remove duplicate photon energies, if any.
    # Duplicate energies can make local derivative estimates unstable,
    # especially near the singular point E' = E.
    #
    # Note:
    # np.unique keeps the first occurrence only because return_index=True.
    # If duplicate energies have different f'' values, this does not average them;
    # it simply keeps the first one.
    E_unique, unique_idx = np.unique(E_grid, return_index=True)
    E_grid = E_unique
    fpp_grid = fpp_grid[unique_idx]

    # At least two points are needed for a numerical integral.
    if len(E_grid) < 2:
        
        raise ValueError("Need at least two energy points for KK integration.")

    # Warn if E_eval is outside the available integration range.
    # np.interp will still return an extrapolated boundary value,
    # but the KK result will be physically unreliable.
    if not (E_grid[0] <= E_eval <= E_grid[-1]):
        
        print(
            f"Warning: E_eval = {E_eval} eV is outside the integration range "
            f"[{E_grid[0]}, {E_grid[-1]}] eV."
        )

    # Define g(E') = E' f''(E').
    # This is the numerator appearing in the KK integrand.
    g_grid = E_grid * fpp_grid

    # Finite integration limits.
    # In the exact KK relation, the integral is from 0 to infinity.
    # Here it is truncated to the available XATOM PACS range.
    E_min = E_grid[0]
    E_max = E_grid[-1]

    # Interpolate f'' to the evaluation energy.
    # This is needed because E_eval may not exactly coincide with a grid point.
    fpp_eval = np.interp(E_eval, E_grid, fpp_grid)

    # g(E_eval) = E_eval * f''(E_eval)
    g_eval = E_eval * fpp_eval

    # Denominator of the KK kernel.
    # This becomes zero when E_grid == E_eval.
    denominator = E_eval**2 - E_grid**2

    # Initialize the subtracted integrand.
    # It will represent:
    #
    #     [g(E') - g(E)] / [E^2 - E'^2]
    #
    # which is finite at E' = E.
    integrand = np.zeros_like(E_grid)

    # Safe division mask.
    # We only divide where the denominator is not zero or extremely tiny.
    # This avoids generating inf/nan from division by zero.
    mask = np.abs(denominator) > 1e-14

    # Numerically evaluate the smooth subtracted part away from the singular point.
    integrand[mask] = (g_grid[mask] - g_eval) / denominator[mask]

    # Handle the possible singular point E' = E_eval.
    # This only happens if E_eval is exactly one of the grid points.
    singular_idx = np.where(~mask)[0]

    if len(singular_idx) > 0:
        idx = singular_idx[0]

        # Estimate g'(E_eval) by finite differences.
        # If the singular point is inside the grid, use a central difference.
        # If it is at the boundary, use one-sided differences.
        if 0 < idx < len(E_grid) - 1:
            g_prime_eval = (
                g_grid[idx + 1] - g_grid[idx - 1]
            ) / (
                E_grid[idx + 1] - E_grid[idx - 1]
            )
        elif idx == 0:
            g_prime_eval = (
                g_grid[idx + 1] - g_grid[idx]
            ) / (
                E_grid[idx + 1] - E_grid[idx]
            )
        else:
            g_prime_eval = (
                g_grid[idx] - g_grid[idx - 1]
            ) / (
                E_grid[idx] - E_grid[idx - 1]
            )

        # Limiting value of the subtracted integrand at E' = E:
        #
        #     lim_{E' -> E} [g(E') - g(E)] / [E^2 - E'^2]
        #     = -g'(E) / (2E)
        #
        # This replaces the singular point with its finite mathematical limit.
        integrand[idx] = -g_prime_eval / (2.0 * E_eval)

    # Numerically integrate the smooth subtracted part.
    # np.trapezoid is valid for non-uniform E_grid because it uses the actual
    # spacing between neighboring energy points.
    numerical_part = np.trapezoid(integrand, E_grid)

    # Analytical principal-value integral of the singular kernel:
    #
    #     P ∫ dE' / (E^2 - E'^2)
    #
    # over the finite range [E_min, E_max].
    #
    # This term accounts for the part that was subtracted from the numerator.
    analytical_pv = (
        1.0 / (2.0 * E_eval)
        * np.log(
            np.abs(
                (E_eval + E_max) * (E_eval - E_min)
                / ((E_eval - E_max) * (E_eval + E_min))
            )
        )
    )

    # Combine the numerical smooth part and the analytical singular-kernel part.
    fp_value = (2.0 / np.pi) * (numerical_part + g_eval * analytical_pv)

    # Return with an extra minus sign to match the negative-imaginary convention
    # used elsewhere in this code.
    return -float(fp_value)





def _cleanup_xatom_files(folder):
    """
    Delete generated XATOM .out and .dat files inside folder.
    """
    folder = Path(folder)
    if not folder.exists():
        return
    for pattern in ("*.out", "*.dat"):
        for file in folder.glob(pattern):
            try:
                file.unlink()
            except FileNotFoundError:
                pass



# def calculate_formfac(config, photonE, temp_storage_folder):
    
#     # Convert XATOM configuration string into an OrderedDict.
#     # Example:
#     #   "3d4,6" becomes:
#     #   "3d-" : 4
#     #   "3d+" : 6
#     #
#     # This is necessary because XATOM's relativistic orbitals are split
#     # into minus and plus components in the output.
#     config_dict = config2dict(config)

#     # Total number of electrons in the current electronic configuration.
#     # This is used as the non-resonant Thomson-like f0 contribution.
#     total_electrons = sum(config_dict.values())

#     # f0 is the normal elastic scattering contribution from all electrons.
#     f0 = total_electrons

#     # fp and fpp accumulate the dispersive and absorptive corrections.
#     #
#     # In this code, the convention is:
#     #
#     #     f = f0 + fp + i fpp
#     #
#     # but absorption/resonance contributes negative fpp.
#     fp, fpp = 0, 0
    
#     # Each photon energy gets its own temporary folder.
#     # This is important when running in parallel, because different workers
#     # should not write into the same XATOM output directory.
#     destination = temp_storage_folder

#     try:
#         # ============================================================
#         # Part 1:
#         # Bound-bound resonant excitation contribution
#         #
#         # This corresponds to discrete intermediate states:
#         #
#         #     i -> n
#         #
#         # where the excited electron remains in a bound orbital.
#         #
#         # This part is evaluated using a Lorentzian resonance formula.
#         # ============================================================

#         # Run XATOM in resonance mode around the current photon energy.
#         # The range photonE ± 50 eV is used to find nearby bound-bound
#         # resonances and transition matrix elements.
#         output_file_res, result_res = run_xatom(
#             config,
#             output_dir=destination,
#             PE_range=f"{photonE-50}-{photonE+50}",
#             is_resonance=True
#         )

#         # Extract all possible bound-bound transitions from the XATOM output.
#         # The resulting table contains:
#         #   from_orb
#         #   to_orb
#         #   E_trans_eV
#         #   dipole_matrix_element
#         bound_bound_transitions = extract_bound_bound_transitions(output_file_res)
        
#         # Loop over all bound-bound transitions found by XATOM.
#         for _, row in bound_bound_transitions.iterrows():

#             # Initial orbital of the excited electron.
#             # Example: "3d+"
#             from_orb = row["from_orb"]

#             # Final bound orbital after excitation.
#             # Example: "5p+"
#             to_orb = row["to_orb"]

#             # Transition energy in eV.
#             E_trans_eV = row["E_trans_eV"]

#             # XATOM transition matrix element is assumed to be in atomic units.
#             # Multiplying by a0 converts it to SI length units in meters.
#             dipole_matrix_element = row["dipole_matrix_element"]*a0

#             # Create a copy of the current configuration so that we can build
#             # the intermediate excited configuration.
#             copy_config_dict = copy.deepcopy(config_dict)

#             # Remove one electron from the initial orbital.
#             copy_config_dict[from_orb] -= 1

#             # Add one electron to the excited bound orbital.
#             copy_config_dict[to_orb] += 1

#             # Convert the modified dictionary back to XATOM config-string format.
#             config_after_transition = dict2config(copy_config_dict)

#             # Run XATOM for the intermediate excited configuration.
#             # This is needed to estimate the lifetime / decay width Gamma_j
#             # of the intermediate state.
#             output_file_excited_state, result_excited_state = run_xatom(
#                 config_after_transition,
#                 output_dir=destination,
#                 PE_range=f"{photonE-1}-{photonE+1}",
#                 is_resonance=True
#             )

#             # Extract total decay rate and lifetime of the intermediate state.
#             decay_info = extract_total_decay_info(output_file_excited_state)

#             # Lifetime in femtoseconds.
#             lifetime_j = decay_info['total_lifetime_fs']

#             # Convert lifetime to decay width in Joules:
#             #
#             #     Gamma = hbar / tau
#             #
#             # lifetime_j is in fs, so multiply by 1e-15 to get seconds.
#             Gamma_j = hbar/(lifetime_j*1e-15)

#             # Real part of Lorentzian resonance contribution.
#             #
#             # This corresponds to the dispersive part f'.
#             #
#             # The numerator contains:
#             #   E - E_trans
#             #
#             # so fp changes sign across resonance.
#             fp_i2j = (2*np.pi*alpha*(photonE*e)**2)/(r_elec*h*c)*dipole_matrix_element**2*(photonE*e-E_trans_eV*e)/((photonE*e-E_trans_eV*e)**2+(Gamma_j/2)**2)

#             # Imaginary part of Lorentzian resonance contribution.
#             #
#             # The minus sign is intentional:
#             # in this code, absorption/resonance corresponds to fpp < 0.
#             fpp_i2j = -(2*np.pi*alpha*(photonE*e)**2)/(r_elec*h*c)*dipole_matrix_element**2*(Gamma_j/2)/((photonE*e-E_trans_eV*e)**2+(Gamma_j/2)**2)

#             # Add this transition's contribution to the total correction.
#             fp += fp_i2j
#             fpp += fpp_i2j

#         # ============================================================
#         # Part 2:
#         # Photoionization contribution
#         #
#         # This corresponds to continuum intermediate states:
#         #
#         #     i -> free electron
#         #
#         # Instead of summing individual continuum states, use the XATOM
#         # photoabsorption cross section and convert it to f''.
#         # Then calculate f' from f'' using the KK relation.
#         # ============================================================

#         # Run XATOM in non-resonant mode over a broad photon-energy range.
#         # The wide range is important for the KK integral because f'(E)
#         # depends on f''(E') over all energies.
#         output_file_nonres, result_nonres = run_xatom(
#             config,
#             output_dir=destination,
#             PE_range=f"0-4000",
#             is_resonance=False
#         )

#         # Parse the full photoabsorption cross section table.
#         PACS_pdframe = extract_pacs_table(output_file_nonres)

#         # Total photoabsorption cross section in Mb.
#         PACS_ionization_total = PACS_pdframe['Total'].to_numpy()

#         # Photon-energy grid corresponding to the PACS values.
#         PhotonEs = PACS_pdframe['P.E.(eV)'].to_numpy()

#         # Convert PACS to f'' using:
#         #
#         #     f'' = E * sigma_abs / 69.9
#         #
#         # but with a negative sign because this code uses the convention
#         # where absorption gives negative imaginary scattering amplitude.
#         fpp_ionization_total = -PACS_ionization_total*PhotonEs/69.9

#         # Calculate the corresponding f' at the requested photon energy
#         # using the KK relation.
#         #
#         # The function internally returns the sign convention compatible
#         # with this negative-fpp convention.
#         fp_ionization_E = kk_fp_single_from_fpp(
#             PhotonEs,
#             fpp_ionization_total,
#             photonE
#         )

#         # Diagnostic printout:
#         # bound-bound fp plus ionization fp.
#         print('f\':', f'{fp}+{fp_ionization_E}')

#         # Diagnostic printout:
#         # bound-bound fpp plus ionization fpp at exactly photonE.
#         #
#         # Note:
#         # This uses exact equality PhotonEs == photonE.
#         # It works only if photonE is exactly present in the XATOM table.
#         print('f\":', f'{fpp}+{fpp_ionization_total[PhotonEs==photonE]}')

#         # Add ionization f' contribution.
#         fp += fp_ionization_E

#         # Add ionization f'' contribution at the requested photon energy.
#         #
#         # Note:
#         # This line returns an array if PhotonEs == photonE matches one entry.
#         # It is later converted to float in the return statement.
#         fpp += fpp_ionization_total[PhotonEs==photonE]

#         # Return:
#         #   f0  : normal electron-count scattering factor
#         #   fp  : real dispersion correction
#         #   fpp : signed imaginary correction
#         return float(f0), float(fp), float(fpp)
    
#     finally:
#         # Always clean up temporary XATOM output files, even if the calculation
#         # fails midway. This prevents parallel runs from leaving many .out/.dat
#         # files behind.
#         _cleanup_xatom_files(destination)

# ======================================================================================================================
# SQLite cache layer
# ======================================================================================================================
# This block contains all database-related helper functions.
#
# Purpose:
#   - Cache photon-energy-independent XATOM results.
#   - Avoid repeating expensive XATOM calls for the same configuration.
#   - Support safe parallel access when many photon energies are calculated
#     at the same time.
#
# Cached quantities:
#   1. Bound-bound transition search results for each original config.
#   2. Excited-state decay lifetime for each after-transition config.
#   3. Non-resonant PE-dependent f'' table for each original config.
#
# Important:
#   Each parallel process should open its own SQLite connection.
#   The database uses WAL mode and a lock table to avoid write conflicts.
# ======================================================================================================================


def config_hash(config):
    """
    Generate a stable hash key for an XATOM configuration string.

    The original config string is still stored in the database for readability,
    but the hash is used as the compact primary key.
    """
    return hashlib.sha1(config.encode("utf-8")).hexdigest()


def _serialize_np_arrays(**arrays):
    """
    Serialize NumPy arrays into a compressed binary blob for SQLite storage.
    """
    buffer = io.BytesIO()
    np.savez_compressed(buffer, **arrays)
    return buffer.getvalue()


def _deserialize_np_arrays(blob):
    """
    Read NumPy arrays back from a SQLite binary blob.

    The arrays are copied before returning so that they no longer depend
    on the internal NpzFile object.
    """
    buffer = io.BytesIO(blob)

    with np.load(buffer) as data:
        arrays = {key: data[key].copy() for key in data.files}

    return arrays


def _open_sqlite_db(db_path):
    """
    Open one SQLite database file with WAL mode enabled.

    Each process should open its own connection.
    """

    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(
        db_path,
        timeout=60.0,
        isolation_level=None,
    )

    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.execute("PRAGMA busy_timeout=60000;")

    conn.execute("""
        CREATE TABLE IF NOT EXISTS cache_locks (
            lock_name TEXT NOT NULL,
            config_hash TEXT NOT NULL,
            owner TEXT NOT NULL,
            created_at REAL NOT NULL,
            PRIMARY KEY (lock_name, config_hash)
        )
    """)

    return conn


def open_bound_cache_db(db_path=BOUND_CACHE_DB_PATH):
    conn = _open_sqlite_db(db_path)

    conn.execute("""
        CREATE TABLE IF NOT EXISTS bound_transition_done (
            config_hash TEXT PRIMARY KEY,
            config TEXT NOT NULL,
            created_at REAL NOT NULL
        )
    """)

    conn.execute("""
        CREATE TABLE IF NOT EXISTS bound_transition_cache (
            config_hash TEXT NOT NULL,
            config TEXT NOT NULL,
            from_orb TEXT NOT NULL,
            to_orb TEXT NOT NULL,
            E_trans_eV REAL NOT NULL,
            dipole_matrix_element REAL NOT NULL,
            config_after_transition TEXT NOT NULL,
            after_transition_lifetime_fs REAL NOT NULL
        )
    """)

    conn.execute("""
        CREATE INDEX IF NOT EXISTS idx_bound_transition_config_hash
        ON bound_transition_cache(config_hash)
    """)

    return conn


def open_decay_cache_db(db_path=DECAY_CACHE_DB_PATH):
    conn = _open_sqlite_db(db_path)

    conn.execute("""
        CREATE TABLE IF NOT EXISTS decay_lifetime_cache (
            config_hash TEXT PRIMARY KEY,
            config TEXT NOT NULL,
            lifetime_fs REAL NOT NULL,
            decay_info_json TEXT,
            created_at REAL NOT NULL
        )
    """)

    return conn


def open_nonres_cache_db(db_path=NONRES_CACHE_DB_PATH):
    conn = _open_sqlite_db(db_path)

    conn.execute("""
        CREATE TABLE IF NOT EXISTS nonres_fpp_cache (
            config_hash TEXT PRIMARY KEY,
            config TEXT NOT NULL,
            data_blob BLOB NOT NULL,
            created_at REAL NOT NULL
        )
    """)

    return conn


def _try_acquire_cache_lock(conn, lock_name, config, timeout=3600, wait=0.2, stale_after=7200):
    """
    Acquire a per-config SQLite lock.

    This avoids the common parallel problem:
        process A checks cache: missing
        process B checks cache: missing
        both run XATOM
        both write

    This function uses an INSERT into a table with a PRIMARY KEY.
    Only one process can insert the same (lock_name, config_hash).

    If a process crashes and leaves a stale lock, the lock is removed after
    stale_after seconds.
    """

    h = config_hash(config)
    owner = f"{os.getpid()}_{time.time()}"
    start = time.time()

    while True:
        now = time.time()

        # Remove stale lock if it is very old.
        conn.execute(
            """
            DELETE FROM cache_locks
            WHERE lock_name = ?
              AND config_hash = ?
              AND created_at < ?
            """,
            (lock_name, h, now - stale_after)
        )

        try:
            conn.execute(
                """
                INSERT INTO cache_locks
                (lock_name, config_hash, owner, created_at)
                VALUES (?, ?, ?, ?)
                """,
                (lock_name, h, owner, now)
            )
            return owner

        except sqlite3.IntegrityError:
            if time.time() - start > timeout:
                raise TimeoutError(
                    f"Timeout waiting for cache lock: "
                    f"lock_name={lock_name}, config_hash={h}"
                )

            time.sleep(wait)


def _release_cache_lock(conn, lock_name, config, owner):
    """
    Release a cache lock.

    The owner check prevents one process from accidentally deleting another
    process's lock.
    """
    h = config_hash(config)

    conn.execute(
        """
        DELETE FROM cache_locks
        WHERE lock_name = ?
          AND config_hash = ?
          AND owner = ?
        """,
        (lock_name, h, owner)
    )


def _load_bound_transitions_from_db(conn, config):
    """
    Load cached bound-bound transitions for one config.

    Returns None if this config has not been marked as completed.
    Returns an empty DataFrame if the config was completed but has no transitions.
    """

    h = config_hash(config)

    done = conn.execute(
        """
        SELECT 1 FROM bound_transition_done
        WHERE config_hash = ?
        """,
        (h,)
    ).fetchone()

    if done is None:
        return None

    rows = conn.execute(
        """
        SELECT
            from_orb,
            to_orb,
            E_trans_eV,
            dipole_matrix_element,
            config_after_transition,
            after_transition_lifetime_fs
        FROM bound_transition_cache
        WHERE config_hash = ?
        """,
        (h,)
    ).fetchall()

    return pd.DataFrame(
        rows,
        columns=[
            "from_orb",
            "to_orb",
            "E_trans_eV",
            "dipole_matrix_element",
            "config_after_transition",
            "after_transition_lifetime_fs",
        ],
    )


def get_decay_lifetime_cached(conn, config_after_transition, output_dir, photonE):
    """
    Case 2 cache:
    Get lifetime of an excited-state config.

    Key:
        config_after_transition

    Cached value:
        lifetime_fs

    If missing, run XATOM once, parse decay information, and save to SQLite.
    """

    h = config_hash(config_after_transition)

    row = conn.execute(
        """
        SELECT lifetime_fs, decay_info_json
        FROM decay_lifetime_cache
        WHERE config_hash = ?
        """,
        (h,)
    ).fetchone()

    if row is not None:
        return float(row[0])

    owner = _try_acquire_cache_lock(
        conn,
        lock_name="decay_lifetime",
        config=config_after_transition,
    )

    try:
        # Check again after acquiring the lock.
        # Another process may have finished while this process was waiting.
        row = conn.execute(
            """
            SELECT lifetime_fs, decay_info_json
            FROM decay_lifetime_cache
            WHERE config_hash = ?
            """,
            (h,)
        ).fetchone()

        if row is not None:
            return float(row[0])

        output_file_excited_state, result_excited_state = run_xatom(
            config_after_transition,
            output_dir=output_dir,
            PE_range=f"{photonE-1}-{photonE+1}",
            is_resonance=True
        )

        decay_info = extract_total_decay_info(output_file_excited_state)
        lifetime_fs = float(decay_info["total_lifetime_fs"])

        conn.execute(
            """
            INSERT OR IGNORE INTO decay_lifetime_cache
            (config_hash, config, lifetime_fs, decay_info_json, created_at)
            VALUES (?, ?, ?, ?, ?)
            """,
            (
                h,
                config_after_transition,
                lifetime_fs,
                json.dumps(decay_info),
                time.time(),
            )
        )

        return lifetime_fs

    finally:
        _release_cache_lock(
            conn,
            lock_name="decay_lifetime",
            config=config_after_transition,
            owner=owner,
        )


def get_bound_transitions_cached(bound_conn, decay_conn, config, output_dir, photonE):
    """
    Case 1 cache:
    Get all possible bound-bound transitions for a given original config.

    Key:
        config

    Cached values:
        config
        from_orb
        to_orb
        E_trans_eV
        dipole_matrix_element
        config_after_transition
        after_transition_lifetime_fs

    The transition search is photon-energy independent because PE_range
    is fixed to 490-770 eV.

    This function uses two independent SQLite connections:
        bound_conn:
            reads/writes the bound-transition cache
        decay_conn:
            reads/writes the excited-state lifetime cache

    Separating these two databases reduces SQLite writer-lock contention
    in parallel SLURM array jobs.
    """

    h = config_hash(config)

    cached_df = _load_bound_transitions_from_db(bound_conn, config)
    if cached_df is not None:
        return cached_df

    owner = _try_acquire_cache_lock(
        bound_conn,
        lock_name="bound_transitions",
        config=config,
    )

    try:
        # Check again after acquiring lock.
        # Another process may have finished this config while this process
        # was waiting for the bound-transition lock.
        cached_df = _load_bound_transitions_from_db(bound_conn, config)
        if cached_df is not None:
            return cached_df

        output_file_res, result_res = run_xatom(
            config,
            output_dir=output_dir,
            PE_range="490-770",
            is_resonance=True
        )

        bound_bound_transitions = extract_bound_bound_transitions(output_file_res)

        config_dict = config2dict(config)
        rows_to_save = []

        for _, row in bound_bound_transitions.iterrows():
            from_orb = row["from_orb"]
            to_orb = row["to_orb"]
            E_trans_eV = float(row["E_trans_eV"])
            dipole_matrix_element = float(row["dipole_matrix_element"])

            copy_config_dict = copy.deepcopy(config_dict)
            copy_config_dict[from_orb] -= 1
            copy_config_dict[to_orb] += 1

            config_after_transition = dict2config(copy_config_dict)

            # The lifetime belongs to the excited-state configuration.
            # It is cached in the independent decay_lifetime_cache.sqlite file.
            after_transition_lifetime_fs = get_decay_lifetime_cached(
                decay_conn,
                config_after_transition=config_after_transition,
                output_dir=output_dir,
                photonE=photonE,
            )

            rows_to_save.append(
                (
                    h,
                    config,
                    from_orb,
                    to_orb,
                    E_trans_eV,
                    dipole_matrix_element,
                    config_after_transition,
                    float(after_transition_lifetime_fs),
                )
            )

        if rows_to_save:
            bound_conn.executemany(
                """
                INSERT INTO bound_transition_cache
                (
                    config_hash,
                    config,
                    from_orb,
                    to_orb,
                    E_trans_eV,
                    dipole_matrix_element,
                    config_after_transition,
                    after_transition_lifetime_fs
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                rows_to_save
            )

        # Mark this config as completed even if it has zero transitions.
        # This avoids repeating a transition search that correctly found no rows.
        bound_conn.execute(
            """
            INSERT OR IGNORE INTO bound_transition_done
            (config_hash, config, created_at)
            VALUES (?, ?, ?)
            """,
            (h, config, time.time())
        )

        return _load_bound_transitions_from_db(bound_conn, config)

    finally:
        _release_cache_lock(
            bound_conn,
            lock_name="bound_transitions",
            config=config,
            owner=owner,
        )

def get_nonres_fpp_cached(conn, config, output_dir):
    """
    Case 3 cache:
    Get PE-dependent non-resonant f'' for one config.

    Key:
        config

    Cached values:
        PhotonEs
        fpp_ionization_total

    The full XATOM .out file is not saved. Only the compact extracted arrays
    are stored in SQLite as a compressed binary blob.
    """

    h = config_hash(config)

    row = conn.execute(
        """
        SELECT data_blob
        FROM nonres_fpp_cache
        WHERE config_hash = ?
        """,
        (h,)
    ).fetchone()

    if row is not None:
        arrays = _deserialize_np_arrays(row[0])
        return arrays["PhotonEs"], arrays["fpp_ionization_total"]

    owner = _try_acquire_cache_lock(
        conn,
        lock_name="nonres_fpp",
        config=config,
    )

    try:
        # Check again after acquiring lock.
        row = conn.execute(
            """
            SELECT data_blob
            FROM nonres_fpp_cache
            WHERE config_hash = ?
            """,
            (h,)
        ).fetchone()

        if row is not None:
            arrays = _deserialize_np_arrays(row[0])
            return arrays["PhotonEs"], arrays["fpp_ionization_total"]

        output_file_nonres, result_nonres = run_xatom(
            config,
            output_dir=output_dir,
            PE_range="0-4000",
            is_resonance=False
        )

        PACS_pdframe = extract_pacs_table(output_file_nonres)

        PACS_ionization_total = PACS_pdframe["Total"].to_numpy()
        PhotonEs = PACS_pdframe["P.E.(eV)"].to_numpy()

        # Signed convention used by your code:
        # absorption gives f'' < 0.
        fpp_ionization_total = -PACS_ionization_total * PhotonEs / 69.9

        data_blob = _serialize_np_arrays(
            PhotonEs=PhotonEs,
            fpp_ionization_total=fpp_ionization_total,
        )

        conn.execute(
            """
            INSERT OR IGNORE INTO nonres_fpp_cache
            (config_hash, config, data_blob, created_at)
            VALUES (?, ?, ?, ?)
            """,
            (h, config, data_blob, time.time())
        )

        return PhotonEs, fpp_ionization_total

    finally:
        _release_cache_lock(
            conn,
            lock_name="nonres_fpp",
            config=config,
            owner=owner,
        )



# ======================================================================================================================
# Main form-factor calculation
# ======================================================================================================================
# This block combines cached XATOM results to calculate:
#
#     f(E) = f0 + f'(E) + i f''(E)
#
# for one electronic configuration and one assigned photon energy.
#
# The expensive photon-energy-independent XATOM calculations are handled
# by the SQLite cache functions above. This function mainly:
#
#   1. Loads or computes bound-bound transition information.
#   2. Adds Lorentzian resonant bound-bound contributions to f' and f''.
#   3. Loads or computes non-resonant ionization f''.
#   4. Calculates ionization f' from f'' using the KK relation.
#   5. Returns f0, f', and f''.
#
# Sign convention:
#   The code uses f = f0 + f' + i f'', with absorption stored as f'' < 0.
# ======================================================================================================================


def calculate_formfac(
    config,
    photonE,
    temp_storage_folder,
    bound_cache_db_path=BOUND_CACHE_DB_PATH,
    decay_cache_db_path=DECAY_CACHE_DB_PATH,
    nonres_cache_db_path=NONRES_CACHE_DB_PATH,
):
    """
    Calculate f0, f', and f'' for one electronic configuration and one photon energy.

    This version uses three separate SQLite cache databases:

        1. bound_transition_cache.sqlite
        2. decay_lifetime_cache.sqlite
        3. nonres_fpp_cache.sqlite

    This reduces write-lock contention compared with storing all cache tables
    inside one SQLite file.
    """

    bound_conn = open_bound_cache_db(bound_cache_db_path)
    decay_conn = open_decay_cache_db(decay_cache_db_path)
    nonres_conn = open_nonres_cache_db(nonres_cache_db_path)

    try:
        config_dict = config2dict(config)
        total_electrons = sum(config_dict.values())

        f0 = total_electrons
        fp, fpp = 0.0, 0.0

        destination = os.path.join(temp_storage_folder, f"{photonE}eV")
        Path(destination).mkdir(parents=True, exist_ok=True)

        try:
            # ============================================================
            # Case 1:
            # Bound-bound transitions, cached by original config.
            # Uses bound_transition_cache.sqlite.
            # The after-transition lifetimes are read/written through
            # decay_lifetime_cache.sqlite.
            # ============================================================

            bound_bound_transitions = get_bound_transitions_cached(
                bound_conn,
                decay_conn,
                config=config,
                output_dir=destination,
                photonE=photonE,
            )

            for _, row in bound_bound_transitions.iterrows():
                E_trans_eV = float(row["E_trans_eV"])
                dipole_matrix_element = float(row["dipole_matrix_element"]) * a0
                lifetime_j = float(row["after_transition_lifetime_fs"])

                Gamma_j = hbar / (lifetime_j * 1e-15)

                E_photon_J = photonE * e
                E_trans_J = E_trans_eV * e

                prefactor = (
                    (2 * np.pi * alpha * E_photon_J**2)
                    / (r_elec * h * c)
                    * dipole_matrix_element**2
                )

                denominator = (
                    (E_photon_J - E_trans_J)**2
                    + (Gamma_j / 2)**2
                )

                fp_i2j = (
                    prefactor
                    * (E_photon_J - E_trans_J)
                    / denominator
                )

                fpp_i2j = (
                    -prefactor
                    * (Gamma_j / 2)
                    / denominator
                )

                fp += fp_i2j
                fpp += fpp_i2j

            # ============================================================
            # Case 3:
            # Non-resonant ionization f'', cached by original config.
            # Uses nonres_fpp_cache.sqlite.
            # ============================================================

            PhotonEs, fpp_ionization_total = get_nonres_fpp_cached(
                nonres_conn,
                config=config,
                output_dir=destination,
            )

            fp_ionization_E = kk_fp_single_from_fpp(
                PhotonEs,
                fpp_ionization_total,
                photonE
            )

            fpp_ionization_E = float(
                np.interp(
                    photonE,
                    PhotonEs,
                    fpp_ionization_total
                )
            )

            print("f':", f"{fp}+{fp_ionization_E}")
            print('f":', f"{fpp}+{fpp_ionization_E}")

            fp += fp_ionization_E
            fpp += fpp_ionization_E

            return float(f0), float(fp), float(fpp)

        finally:
            _cleanup_xatom_files(destination)

    finally:
        bound_conn.close()
        decay_conn.close()
        nonres_conn.close()