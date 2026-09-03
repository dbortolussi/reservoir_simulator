import os
import sys
import numpy as np

lib_path = r"C:\Users\Daniel\Documents\C++\reservoir_project\out\build\x64-release"
sys.path.insert(0, lib_path)

import pyressim

def make_simulator(K, mu_w=1.0e-3, mu_o=5.0e-3, phi=0.2, Lx=100.0, Ly=100.0, p_inj=3.0e7, p_prod=1.0e7):
    nx, ny = K.shape
    sim = pyressim.Simulator(nx, ny, Lx, Ly)
    sim.set_permeability(np.ascontiguousarray(K, dtype=float))
    sim.set_fluid(mu_w=mu_w, mu_o=mu_o)
    sim.set_relperm(Swc=0.2, Sor=0.2, n=2.0)
    sim.set_porosity(phi)
    sim.add_well(0,      0,      "injector", p_inj)
    sim.add_well(nx - 1, ny - 1, "producer", p_prod)
    return sim

def summarize(result, label=""):
    it = np.array(result.newton_iters, dtype=float)
    mean_it = it.mean() if it.size else 0.0
    max_it = int(it.max()) if it.size else 0
    print(f"[{label}] steps={len(result.dt_history)}  "
    f"rollbacks={int(np.sum(result.rollbacks))}  "
    f"mean_newton={mean_it:.2f}  max_newton={max_it}  "
    f"breakthrough={result.breakthrough_time:.3e} s")

