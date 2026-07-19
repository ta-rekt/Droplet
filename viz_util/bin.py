"""
load_proximal_promoter.py

Configuration + loader for the proximal_promoter_transc_d_50 experiment.
Merges runs from 20260407 and 20260413 (16 runs each) at each concentration.

Paste this into a notebook cell and run it. Produces:
    on_times, off_times, bound_times, unbound_times,
    pattern_distances, transcription_activity, binding_activity,
    local_concentration   -- each a list of length k (parameters),
                             where each entry is a list of m numpy arrays.
    close_mask            -- boolean version of pattern_distances < 5,
                             same shape.
    params                -- list of parameter values (concentrations).
"""

import importlib
import viz_util.load_transcription as lt
importlib.reload(lt)

# Concentration values shared across both date folders
_concs = [
    0.003, 0.004, 0.005, 0.006, 0.007, 0.008,
    0.009, 0.010, 0.011, 0.012, 0.013, 0.014,
    0.015, 0.016, 0.017, 0.018,
]
_dates = ["data/20260423"]
_exp   = "proximal_promoter_transc_d_50_c_{:.3f}"
_runs  = range(1, 17)   # output_cpu1.json ... output_cpu16.json

groups = lt.make_groups([
    (c, [("{}/{}".format(d, _exp.format(c)), _runs) for d in _dates])
    for c in _concs
])

(on_times, off_times, bound_times, unbound_times,
pattern_distances, transcription_activity, binding_activity,
local_concentration, run_ids) = lt.load(groups, n1=0, n=0)

close_mask = lt.make_mask(pattern_distances, threshold=5.0)

params = _concs

print("Loaded {} parameter values, {} runs each.".format(
    len(params), len(on_times[0])))


####### plotting 

# average across independent runs

import numpy as np
import matplotlib.pyplot as plt

nruns = 16
nconc = 16

tf_conc = np.linspace(0.003, 0.018, 16)
mean_transcription = np.zeros(nconc)
mean_local_concentration = np.zeros(nconc)
mean_pattern_distances = np.zeros(nconc)
stuck_probability = np.zeros(nconc)

stdev_transcription = np.zeros(nconc)
stdev_local_concentration = np.zeros(nconc)
stdev_pattern_distances = np.zeros(nconc)

for cidx in range(nconc):
    mean_transcription[cidx] = np.mean(
        [np.mean(transcription_activity[cidx][run]) for run in range(nruns)]
    )
    mean_local_concentration[cidx] = np.mean(
        [np.mean(local_concentration[cidx][run]) for run in range(nruns)]
    )
    mean_pattern_distances[cidx] = np.mean(
        [np.mean(pattern_distances[cidx][run]) for run in range(nruns)]
    )
    stuck_probability[cidx] = np.mean(
        [np.mean(pattern_distances[cidx][run] < 4) for run in range(nruns)]
    )

    stdev_transcription[cidx] = np.std(
        [np.mean(transcription_activity[cidx][run]) for run in range(nruns)]
    )
    stdev_local_concentration[cidx] = np.std(
        [np.mean(local_concentration[cidx][run]) for run in range(nruns)]
    )
    stdev_pattern_distances[cidx] = np.std(
        [np.mean(pattern_distances[cidx][run]) for run in range(nruns)]
    )

fig, ax = plt.subplots(3,1, figsize=(5,8))

ax[0].plot(tf_conc, mean_transcription)
for cidx in range(nconc):
    ax[0].scatter(np.ones(nruns)* tf_conc[cidx], 
        [np.mean(transcription_activity[cidx][run]) for run in range(nruns)],
        color="C0", alpha=0.2)

ax[0].set_ylabel("Transcription", fontsize=12)

ax[1].plot(tf_conc, mean_local_concentration, color="C2")
ax[1].plot(tf_conc, tf_conc, color="C2", linestyle="dotted", label="$\phi_{TF}$")
for cidx in range(nconc):
    ax[1].scatter(np.ones(nruns)* tf_conc[cidx], 
        [np.mean(local_concentration[cidx][run]) for run in range(nruns)],
        color="C2", alpha=0.2)

ax[1].set_ylabel("Local concentration", fontsize=12)
ax[1].set_yscale('log')
ax[1].axhline(y=0.05, linestyle="dashed", color="black", alpha=0.3, label="$K_d$")
ax[1].legend()


ax[2].plot(tf_conc, mean_pattern_distances, color="C1")
for cidx in range(nconc):
    ax[2].scatter(np.ones(nruns)* tf_conc[cidx], 
        [np.mean(pattern_distances[cidx][run]) for run in range(nruns)],
        color="C1", alpha=0.2)

ax[2].set_ylabel("Pattern distances", fontsize=12)
ax[2].set_xlabel("Concentration $\phi_{TF}$", fontsize=12)