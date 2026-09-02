import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
import numpy as np
from permeability_fields import lognormal_field
from simcommon import make_simulator
import matplotlib.pyplot as plt

def main():
    nx = ny = 48

    t_total = 1.5e7

    sigmas = np.linspace(-4, 2, num=100)

    rerun = True

    if rerun:

        results = []

        for sig in sigmas:

            print(f"Running sigma {sig}")
    
            K = lognormal_field(nx, ny, sigma=10**sig)

            sim = make_simulator(K)
            #try:
            #    res_newton = sim.run("newton", t_end=t_total)
            #except:
            #    print("failed")
            #    continue
            res_impes = sim.run("impes", t_end=t_total)

            result=res_impes
            nx, ny = result.nx, result.ny
            Sw = np.array(result.Sw_final).reshape(nx, ny)
            plt.imsave(rf'pictures\res_{str(sig).replace("-", "minus")}.png', Sw)

            #mae = np.median(np.abs(res_newton.Sw_final - res_impes.Sw_final))

            results.append((sig, sum(res_impes.newton_iters)))
        
        np.save('results.npy', results)
    
    else:
        results = np.load('results.npy')

    results = np.array(results)
    plt.semilogy(results[:, 0], results[:, 1], '.')
    plt.show()


if __name__ == "__main__":
    main()