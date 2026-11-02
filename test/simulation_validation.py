import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
import numpy as np
from permeability_fields import log_correlated_gaussian
import pyressim
import matplotlib.pyplot as plt


def test_1d_known_solution(solver: str):

    # assumes corey-brooks relation

    total_length = 100
    mu_w=1.0e-3
    mu_o=5.0e-3
    phi=0.2
    perm = 1e-13
    injector_p = 3.0e7
    producer_p = 1.0e7
    Swc=0.2
    Sor=0.2
    n=2.0

    # Known solution uses Buckley-Leverett for saturation front halfway 
    # across the length of the interval

    visc_ratio = mu_o / mu_w
    mobile_sat_range = 1 - Swc - Sor
    norm_front_sat = 1 / np.sqrt(1 + visc_ratio)
    saturation_range = np.linspace(1.0, norm_front_sat, num=10001)
    inverse_mobilities = 1 / (saturation_range**2/mu_w + (1-saturation_range)**2/mu_o)
    denom = visc_ratio * saturation_range**2 + (1-saturation_range)**2
    df_dSw = visc_ratio * 2 * saturation_range * (1-saturation_range) / (mobile_sat_range * denom**2)
    front_pos_to_inj_water_dist = (1 + np.sqrt(1+visc_ratio)) / (2 * mobile_sat_range)
    x_fraction = df_dSw / front_pos_to_inj_water_dist
    average_inverse_mobility = np.sum(0.5 * (inverse_mobilities[:-1] + inverse_mobilities[1:]) * np.diff(x_fraction))
    
    x_f = total_length / 2
    a = x_f / front_pos_to_inj_water_dist
    front_transit_time = phi / (perm * (injector_p - producer_p)) * (total_length*mu_o*a + (front_pos_to_inj_water_dist*(average_inverse_mobility-mu_o)/2) * a**2)

    print(f"front transit time {front_transit_time}")

    nx = 200
    ny = 1
    Lx=Ly=total_length/nx

    K = perm * np.ones((nx, ny))
    

    sim = pyressim.Simulator(nx, ny, Lx, Ly)
    sim.set_permeability(np.ascontiguousarray(K, dtype=float))
    sim.set_fluid(mu_w=mu_w, mu_o=mu_o)
    sim.set_relperm(Swc=Swc, Sor=Sor, n=n)
    sim.set_porosity(phi)
    sim.add_well(0, 0, "injector", injector_p)
    sim.add_well(nx-1, 0, "producer", producer_p)

    res = sim.run(solver, t_end=front_transit_time)
    result=res

    nx, ny = result.nx, result.ny
    Sw = np.array(result.Sw_final).reshape(nx, ny)
    plt.plot(np.linspace(0, total_length, num=nx), Sw)
    plt.show()


test_1d_known_solution("newton")
