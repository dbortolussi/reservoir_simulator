import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
import numpy as np
from permeability_fields import log_correlated_gaussian
from simcommon import make_simulator
import matplotlib.pyplot as plt

def main():
    nx = ny = 48

    t_total = 1.5e7

    K = log_correlated_gaussian(nx)

    sim = make_simulator(K)

    print("made simulator")

    res_impes = sim.run("impes", t_end=t_total)

    result=res_impes
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