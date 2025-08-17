import numpy as np

def get_gas_properties(gas: str):
    """
    Returns n2 in m^2/W at 1 atm, and approximate gas refractive index.
    """
    gas = gas.lower()
    gas_data = {
        'argon': {'n2': 1.4e-23, 'n': 1.000281},
        'neon':  {'n2': 0.74e-23, 'n': 1.000067},
        'air':   {'n2': 3.2e-23, 'n': 1.000293},
        'krypton': {'n2': 2.4e-23, 'n': 1.000427},
        'xenon': {'n2': 4.2e-23, 'n': 1.000702},
    }
    if gas not in gas_data:
        raise ValueError(f"Gas '{gas}' not found. Available: {list(gas_data.keys())}")
    return gas_data[gas]['n2'], gas_data[gas]['n']

def calculate_gamma_beta2(
    gas='argon', 
    pressure_bar=1.0, 
    core_radius_um=15.0, 
    wavelength_nm=800.0
):
    """
    Calculate nonlinear parameter γ (W⁻¹·m⁻¹) and GVD β₂ (s²/m)
    """
    # Physical constants
    c = 3e8  # Speed of light (m/s)
    pi = np.pi
    u01 = 2.405  # First zero of Bessel function (HE11 mode)
    
    # Convert units
    wavelength_m = wavelength_nm * 1e-9
    omega = 2 * pi * c / wavelength_m
    core_radius_m = core_radius_um * 1e-6
    A_eff = pi * core_radius_m**2  # Approximation

    # Gas properties
    n2_1atm, n = get_gas_properties(gas)
    n2 = n2_1atm * pressure_bar  # Scale with pressure

    # Nonlinear coefficient γ
    gamma = (2 * pi * n2) / (wavelength_m * A_eff)

    # GVD β₂ (only waveguide dispersion, assuming negligible material dispersion)
    beta2 = - (u01**2) / (2 * pi * c * core_radius_m**2) * (1 / omega**2)

    return gamma, beta2

gamma, beta2 = calculate_gamma_beta2(
    gas='argon', 
    pressure_bar=1.0, 
    core_radius_um=15.0, 
    wavelength_nm=800.0
)

print(f"γ = {gamma:.2e} W⁻¹·m⁻¹")
print(f"β₂ = {beta2:.2e} s²/m")