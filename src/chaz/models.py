import numpy as np


def r_max(v_max: np.ndarray, latitude: np.ndarray) -> np.ndarray:
    """
    Tropical cyclone radius to maximum wind as a function of maximum sustained
    wind speed and latitude.

    Log-linear fit to aircraft observation data.

    Source: https://doi.org/10.1175/MWR2831.1, equation 12.1.

    N.B. The source does not include taking the absolute value of the latitude,
    but I cannot see how it could be correct for both hemispheres otherwise.
    """
    return 46.29 * np.exp(-0.0153 * v_max + 0.0166 * np.abs(latitude))
