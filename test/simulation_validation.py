import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "python"))
import numpy as np
from permeability_fields import log_correlated_gaussian
import pyressim
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def curve_compare(computed_curve: np.array, exact_curve: np.array):

    c1 = computed_curve
    c2 = exact_curve

    if c1.shape[1] != 2 or c2.shape[1] != 2:
        raise ValueError(f"wrong shapes {c1.shape} and {c2.shape}") 

    x1 = c1[:,0]
    x2 = c2[:,0]
    y1 = c1[:,1]
    y2 = c2[:,1]

    x_common = np.linspace(min(min(x1), min(x2)), max(max(x1), max(x2)), 1000)

    f1 = interp1d(x1, y1, kind='linear', fill_value="extrapolate")
    f2 = interp1d(x2, y2, kind='linear', fill_value="extrapolate")

    y1_interp = f1(x_common)
    y2_interp = f2(x_common)

    area_between = np.trapz(np.abs(y1_interp - y2_interp), x_common)
    area_total =  np.trapz(np.abs(y2_interp), x_common)

    return area_between / area_total



def test_1d_known_solution(solver: str):

    num_samples = 500

    # assumes corey-brooks relation
    # returns relative area between computed and exact saturation curves

    total_length = 100
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

    # Known solution uses Buckley-Leverett for saturation front halfway 
    # across the length of the interval

    visc_ratio = mu_o / mu_w
    mobile_sat_range = 1 - Swc - Sor
    norm_front_sat = 1 / np.sqrt(1 + visc_ratio)
    saturation_range = np.linspace(1.0, norm_front_sat, num=num_samples)
    inverse_mobilities = 1 / (saturation_range**2/mu_w + (1-saturation_range)**2/mu_o)
    denom = visc_ratio * saturation_range**2 + (1-saturation_range)**2
    df_dSw = visc_ratio * 2 * saturation_range * (1-saturation_range) / (mobile_sat_range * denom**2)
    front_pos_to_inj_water_dist = (1 + np.sqrt(1+visc_ratio)) / (2 * mobile_sat_range)
    x_fraction = df_dSw / front_pos_to_inj_water_dist
    average_inverse_mobility = np.sum(0.5 * (inverse_mobilities[:-1] + inverse_mobilities[1:]) * np.diff(x_fraction))
    
    x_f = total_length * front_point_frac
    a = x_f / front_pos_to_inj_water_dist
    front_transit_time = phi / (perm * (injector_p - producer_p)) * (total_length*mu_o*a + (front_pos_to_inj_water_dist*(average_inverse_mobility-mu_o)/2) * a**2)
    x_exact = np.r_[x_f * x_fraction, x_f, total_length]
    Sw_exact = np.r_[Swc + mobile_sat_range * saturation_range, Swc, Swc]

    print(f"front transit time {front_transit_time}")

    nx = num_samples
    ny = 1
    Lx=total_length
    Ly=total_length / nx

    K = perm * np.ones((nx, ny))
    

    sim = pyressim.Simulator(nx, ny, Lx, Ly)
    sim.set_permeability(np.ascontiguousarray(K, dtype=float))
    sim.set_fluid(mu_w=mu_w, mu_o=mu_o)
    sim.set_relperm(Swc=Swc, Sor=Sor, n=n)
    sim.set_porosity(phi)
    sim.add_well(0, 0, "injector", injector_p)
    sim.add_well(nx-1, 0, "producer", producer_p)
    sim.dt_init = 3000
    res = sim.run(solver, t_end=front_transit_time)
    result=res

    Sw = np.array(result.Sw_final)
    x_vals = np.linspace(0, Lx, nx) + Lx/nx/2

    plt.plot(x_vals, Sw, '.', label="Computed")
    plt.plot(x_exact, Sw_exact, label="Buckley–Leverett")
    plt.show()

    
    return curve_compare(np.vstack((x_vals, Sw)).T, np.vstack((x_exact, Sw_exact)).T)

solvers = ["newton", "impes"]

for solver in solvers:
    err = test_1d_known_solution(solver)
    print(f"Solver {solver} ran with error {100*err}%")
