import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
import numpy as np
import pyressim
from permeability_fields import log_correlated_gaussian
import matplotlib.pyplot as plt

def main():

    nx = ny = 72

    t_total = 2e7

    Lx=Ly=100
    mu_w=1.0e-3
    mu_o=5.0e-3
    phi=0.5
    injector_p = 3.0e8
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

    
    K = log_correlated_gaussian(nx, sigma=1e-1, seed=int(nx*t_total))

    sim.set_permeability(K)
    # sim.dt_init = 10000
    result = sim.run("newton", t_end=t_total)
    Sw = result.Sw_final

    fig, (ax1, ax2) = plt.subplots(1, 2)

    ax2 .imshow(Sw.reshape((nx, ny)))
    plt.colorbar(ax2.images[0], shrink=0.5)
    ax1.imshow(np.log10(K))
    plt.colorbar(ax1.images[0], shrink=0.5)

    plt.tight_layout()
    plt.show()



if __name__ == "__main__":
    main()