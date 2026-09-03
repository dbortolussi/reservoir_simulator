import numpy as np


def lognormal_field(nx, ny, mean_k=1e-13, sigma=1.0, seed=0):
    """Mildly heterogeneous, spatially-correlated log-normal permeability."""
    rng = np.random.default_rng(seed)
    field = rng.normal(0.0, sigma, size=(nx, ny))
    for _ in range(4):                       # smooth for spatial correlation
        field = 0.25 * (np.roll(field, 1, 0) + np.roll(field, -1, 0)
                        + np.roll(field, 1, 1) + np.roll(field, -1, 1))
    return mean_k * np.exp(field)


def layered_field(nx, ny, k_hi=1e-12, k_lo=1e-15, n_layers=5):
    """Alternating high/low permeability horizontal layers."""
    K = np.empty((nx, ny))
    edges = np.linspace(0, nx, n_layers + 1).astype(int)
    for L, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        K[a:b, :] = k_hi if L % 2 == 0 else k_lo
    return K


def channelized_field(nx, ny, k_bg=1e-15, k_ch=1e-12, seed=1):
    """A meandering high-permeability channel in a tight background."""
    rng = np.random.default_rng(seed)
    K = np.full((nx, ny), k_bg)
    j = ny // 2
    for i in range(nx):
        j = int(np.clip(j + rng.integers(-1, 2), 1, ny - 2))
        K[i, j - 1:j + 2] = k_ch
    return K


def killer_case(nx, ny, k_bg=1e-13, k_barrier=1e-18):
    """High-contrast barrier with a single narrow gap: designed to stress
    the Newton solver via localized, near-discontinuous pressure gradients."""
    K = np.full((nx, ny), k_bg)
    wall = nx // 2
    K[wall, :] = k_barrier                   # near-impermeable vertical wall
    gap = ny // 2
    K[wall, gap - 1:gap + 2] = k_bg          # the only path through the wall
    return K