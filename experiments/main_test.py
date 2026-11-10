import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
import numpy as np
import pyressim
from permeability_fields import log_correlated_gaussian
import matplotlib.pyplot as plt
import time
from scipy.ndimage import label

def timer(sim, type, total_time):
    t1 = time.time()
    sim.run(type, t_end=total_time)
    return time.time() - t1

def channelization_metric(K):

    unraveled = np.ravel(K)
    n = len(unraveled)

    def spanned(matrix):
        rows = np.prod(matrix, axis=0)
        cols = np.prod(matrix, axis=1)
        return any(rows) or any(cols)
        
    for ind, val in enumerate(sorted(unraveled, reverse=True)):
        if spanned(K >= val):
            return (n - ind - 1) / n
    else:
        raise IndexError("how did i get here")

def spanning_fraction(k):
    order = np.argsort(k, axis=None)[::-1]
    lo, hi = 0, k.size
    while hi - lo > 1:
        mid = (lo + hi) // 2
        mask = np.zeros(k.shape, dtype=bool)
        mask.flat[order[:mid]] = True
        components, _ = label(mask)
        common = np.intersect1d(components[0], components[-1])
        if np.any(common > 0):
            hi = mid
        else:
            lo = mid
    return hi / k.size

def main():

    nx = ny = 36

    t_total = 2e7

    Lx=Ly=100
    mu_w=1.0e-3
    mu_o=5.0e-3
    phi=0.2
    perm = 1e-13
    injector_p = 3.0e7
    producer_p = 1.0e7
    Swc=0.2
    Sor=0.2
    n=2.0

    sim = pyressim.Simulator(nx, ny, Lx, Ly)
    sim.set_fluid(mu_w=mu_w, mu_o=mu_o)
    sim.set_relperm(Swc=Swc, Sor=Sor, n=n)
    sim.set_porosity(phi)
    sim.add_well(0, 0, "injector", injector_p)
    sim.add_well(nx-1, ny-1, "producer", producer_p)
    #sim.dt_init = 300

    maxes = []
    water_vols = []
    pvis = []
    c_arr = []
    num_samples = 500
    for seed in range(num_samples):

        K = log_correlated_gaussian(nx, seed=seed, sigma=1e-1, scale=1e-13)

        # plt.imshow(np.log(K))
        # plt.show()

        c = spanning_fraction(K)
        c_arr.append(c)
        # continue

        print(seed, c)
        sim.set_permeability(K)

        t_small = 900
        sim.dt_init = t_small
        result_small = sim.run("newton", t_end=t_total)
        Sw_ref = result_small.Sw_final

        t_l = 1000
        t_h = 4500000
        max_steps = 20
        t = 1000
        while (t := t_l + (t_h - t_l) / 2) and (t_h - t_l > 500) and (max_steps := max_steps - 1):

            print(f"trying t={t} t_h={t_h} t_l={t_l}")
            try:
                sim.dt_init = sim.dt_min = t
                result_newton = sim.run("newton", t_end=t_total)
                Sw = result_newton.Sw_final
                error = np.mean(np.abs(Sw - Sw_ref))
                if np.isfinite(error) and error < .1:
                    t_l = t
                    continue
            except RuntimeError as exc:
                if "Timestep fell below minimum" not in str(exc):
                    raise
            t_h = t
        
        water_vol = result_small.injected_water
        pvi = water_vol / (phi * Lx * Ly)
        maxes.append(t)
        water_vols.append(water_vol)
        pvis.append(pvi)

    np.save('c_arr.npy', c_arr)

    maxes = np.load('maxes.npy')
    pvis = np.load('pvis.npy')
    slope, intercept = np.polyfit(
    np.log10(maxes), np.log10(pvis), 1
)
    print(f"PVI ≈ {10**intercept:.4g} * dt**({slope:.3f})")

    plt.plot(c_arr, pvis * maxes / t_total, '.')

    # np.save('maxes.npy', maxes)
    # np.save('pvis.npy', pvis)

    plt.figure()
    plt.loglog(maxes, pvis, '.')

    plt.show()


if __name__ == "__main__":
    main()