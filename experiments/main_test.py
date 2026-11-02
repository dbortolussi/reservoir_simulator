import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
import numpy as np
import pyressim
from permeability_fields import log_correlated_gaussian
import matplotlib.pyplot as plt

def main():
    nx = ny = 48

    t_total = 6e7

    Lx=Ly=100
    front_point_frac = 0.75
    mu_w=1.0e-3
    mu_o=5.0e-3
    phi=0.2
    perm = 1e-13
    injector_p = 3.0e7
    producer_p = 1.0e7
    Swc=0.2
    Sor=0.2
    n=2.0

    K = log_correlated_gaussian(nx)

    sim = pyressim.Simulator(nx, ny, Lx, Ly)
    sim.set_permeability(np.ascontiguousarray(K, dtype=float))
    sim.set_fluid(mu_w=mu_w, mu_o=mu_o)
    sim.set_relperm(Swc=Swc, Sor=Sor, n=n)
    sim.set_porosity(phi)
    sim.add_well(0, 0, "injector", injector_p)
    sim.add_well(nx-1, 0, "producer", producer_p)
    result = sim.run("newton", t_end=t_total)

    nx, ny = result.nx, result.ny
    Sw = np.array(result.Sw_final).reshape(nx, ny)
    plt.imshow(Sw)
    plt.show()
    #plt.imsave(rf'pictures\res_{str(sig).replace("-", "minus")}.png', Sw)

    #results = np.array(results)
    #plt.semilogy(results[:, 0], results[:, 1], '.')
    #plt.show()


if __name__ == "__main__":
    main()