XATOM Form Factor Calculator

Python workflow for calculating complex atomic scattering factors from XATOM outputs:

[
f(E) = f_0 + f’(E) + i f’’(E)
]

The code was developed for large-scale iodine electronic-configuration calculations, but the structure can be adapted to other atomic systems supported by XATOM.

Features

* Runs XATOM for electronic configurations.
* Parses bound-bound resonant transitions.
* Extracts excited-state decay lifetimes.
* Converts photoabsorption cross sections to (f’’).
* Calculates (f’) from (f’’) using a Kramers-Kronig principal-value integral.
* Uses SQLite cache files to avoid repeated photon-energy-independent XATOM calculations.
* Supports SLURM array jobs for parallel photon-energy scans.

Sign convention

The code uses

[
f(E) = f_0 + f’(E) + i f’’(E)
]

with absorption stored as negative (f’’):

[
f’’ < 0
]

The photoionization contribution is therefore calculated as

fpp_ionization_total = -PACS_ionization_total * PhotonEs / 69.9

where the photoabsorption cross section is in Mb and photon energy is in eV.

Main files

calculate_formfac_xatom.py
Parallel_Computing/
    parallel_simulation.py
    photon_energies.txt
    run_array.slurm

calculate_formfac_xatom.py contains the main XATOM interface, parsers, SQLite cache functions, and form-factor calculation.

parallel_simulation.py runs all configurations for one photon energy.

run_array.slurm submits photon-energy calculations as a SLURM array job.

SQLite cache

The workflow caches three types of photon-energy-independent XATOM results:

Cache/
    bound_transition_cache.sqlite
    decay_lifetime_cache.sqlite
    nonres_fpp_cache.sqlite

These cache files store bound-bound transitions, excited-state lifetimes, and non-resonant (f’’(E)) tables. This prevents repeated XATOM runs for configurations that have already been evaluated.

Usage

Generate photon energies:

seq 540 720 > Parallel_Computing/photon_energies.txt

Run one photon energy manually:

python Parallel_Computing/parallel_simulation.py --photonE 680

Submit the SLURM array:

cd Parallel_Computing
sbatch run_array.slurm

Recommended first test:

#SBATCH --array=1-181%10

This runs 181 photon-energy tasks in total, with at most 10 running simultaneously.

Important setup note

The function run_xatom() contains the path to the local XATOM executable. This path must be changed for a different machine or installation.

The Conda/module environment should be activated before running the Python calculation. In the current SLURM workflow, this is handled inside run_array.slurm.

Output

For each photon energy, results are saved as tab-separated text:

ConfigFiles/<photonE>eV/formfac_<photonE>eV.txt

with columns:

config    f0    fp    fpp

Errors are saved separately:

ConfigFiles/<photonE>eV/formfac_<photonE>eV_errors.txt

This makes interrupted calculations restartable, because already-finished configurations can be skipped on rerun.

Notes

This repository is mainly a research workflow for XATOM-based resonant form-factor calculations. Large generated files, temporary XATOM outputs, SQLite caches, and result folders should usually be excluded from Git tracking.