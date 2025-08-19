import numpy as np

def refractive_index_fused_silica(lambda_nm):
    """
    Calculate the refractive index of fused silica (SiO2) using the Sellmeier equation.
    Input wavelength should be in nanometers.
    """
    # Convert nm to micrometers
    lam_um = lambda_nm / 1000.0

    # Sellmeier coefficients for fused silica
    B1 = 0.6961663
    B2 = 0.4079426
    B3 = 0.8974794
    C1 = 0.0684043**2
    C2 = 0.1162414**2
    C3 = 9.896161**2

    lam2 = lam_um**2

    n_squared = 1 + (B1 * lam2) / (lam2 - C1) + (B2 * lam2) / (lam2 - C2) + (B3 * lam2) / (lam2 - C3)
    n = np.sqrt(n_squared)
    return n

# Example: refractive index at 800 nm
lambda_nm = 800
n_800 = refractive_index_fused_silica(lambda_nm)
print(f"Refractive index of fused silica at {lambda_nm} nm = {n_800:.6f}")
