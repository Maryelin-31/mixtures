"mixture fluids, N2 + O2"

" Specific Heat of nitrogen Gas  from: https://www.engineeringtoolbox.com/nitrogen-d_977.html"

from zmlx.fluid.conf.gas_density.O2_density import *
from zmlx.fluid.conf.gas_viscosity.O2_viscosity import *
from zmlx.fluid.conf.gas_density.n2 import *
from zmlx.fluid.conf.gas_viscosity.N2_viscosity import *
from zml import Interp2, Seepage

import numpy as np
import matplotlib.pyplot as plt


def density_mix (dens, fractions):
    """
    Calculates the density of a gas mixture using a weighted average of the densities of the individual components.

     Parameters:
     densities (list): List of densities of the individual components (in kg/m^3 or any other consistent unit).
     mole_fractions (list): List of mole fractions of the components.
    
     Returns:
     float: Density of the mixture (in the same unit as the individual densities).
    """
    if len(dens) != len(fractions):
        raise ValueError("Lists of densities and mole fractions must be the same length")
    
    # Calculate the density of the mixture
    den_mix = sum(f * d for f, d in zip(fractions, dens))
    
    return den_mix


def viscosity_mix(vis, fractions, molecular_mass):
    """
     Calculate the viscosity of a gas mixture using Wilke's Law.
     https://doi.org/10.1063/1.1747673

     Parameters:
     viscosities (list): List of viscosities of the individual components (in Pa*s).
     mole_fractions (list): List of mole fractions of the components.
     molecular_masses (list): List of molecular masses of the components (in Kg/mol).
    
     Returns:
     float: Viscosity of the mixture (in Pa*s).
    """
    n = len(vis)
    phi = np.zeros((n, n))
   
    # Calculate the interaction factor phi_ij
    for i in range(n):
        for j in range(n):
            if i != j:
               phi[i, j] = (1 + (vis[i] / vis[j])**0.5 * (molecular_mass[j] / molecular_mass[i])**0.25)**2 / (8 * (1 + molecular_mass[i] / molecular_mass[j]))**0.5
            else:
               phi[i, j] = 1
   
   # Calculate the viscosity of the mixture
    vis_mix = 0
    for i in range(n):
       suma_phi = 0
       for j in range(n):
           suma_phi += fractions[j] * phi[i, j]
       vis_mix += fractions[i] * vis[i] / suma_phi             
    return vis_mix


def specific_heat_mix(spec_heats, fractions):
    """
     Calculate the specific heat of a mixture of gases.

     Parameters:
     specific_heats (list): List of specific heats of the individual components (in J/(kg·K)).
     mole_fractions (list): List of mole fractions of the components.
    
     Returns:
     float: Specific heat of the mixture (in J/(kg·K)).
    """
    
    if len(spec_heats) != len(fractions):
        raise ValueError("Las listas de calores específicos y fracciones molares deben tener la misma longitud")
    
    # Calculate the specific heat of the mixture
    sh_mix = sum(f * sh for f, sh in zip(fractions, spec_heats))
    
    return sh_mix

def create_flu(tmin=280, tmax=700, pmin=1.0e6, pmax=20.0e6):
    
    assert 250 < tmin < tmax < 1000
    assert 0.01e6 < pmin < pmax < 30.0e6
    
    def mix_den(P, T):
        density = density_mix(dens=[den_O2(P, T), den_N2(P, T)], fractions=[0.1, 0.9])
        return density
    
    def get_density(P, T):
        return mix_den(P, T)
    
    def create_density():
        den = Interp2()
        den.create(pmin, 1e6, pmax, 280, 10, 700, get_density)
        return den
    
    def mix_vis(P, T):
        viscosity = viscosity_mix(vis=[vis_o2(P, T), vis_N2(P, T)], fractions=[0.1, 0.9], molecular_mass=[0.031999, 0.028013])
        return viscosity 
    
    def get_viscosity(P, T):
        return mix_vis(P, T)
    
    def create_viscosity():
        vis = Interp2()
        vis.create(pmin, 1e6, pmax, tmin, 10, tmax, get_viscosity)
        return vis
    
    specific_heat =  specific_heat_mix(spec_heats=[1090, 1167], fractions=[0.1, 0.9])   

    return Seepage.FluDef(den=create_density(), vis=create_viscosity(), specific_heat=specific_heat)

if __name__ == '__main__':
    flu = create_flu()
    
    # Create pressure and temperature meshes
    P_values = np.linspace(1.0e6, 20.0e6, 100)
    T_values = np.linspace(280, 700, 100)
    P_grid, T_grid = np.meshgrid(P_values, T_values)

    # Calculate density and viscosity as a function of P and T
    density_grid = np.zeros_like(P_grid)
    viscosity_grid = np.zeros_like(P_grid)

    for i in range(P_grid.shape[0]):
        for j in range(P_grid.shape[1]):
            p = P_grid[i, j]
            t = T_grid[i, j]
            density_grid[i, j] = density_mix(dens=[den_O2(p, t), den_N2(p, t)], fractions=[0.1, 0.9])
            viscosity_grid[i, j] = viscosity_mix(vis=[vis_o2(p, t), vis_N2(p, t)], fractions=[0.1, 0.9], molecular_mass=[0.031999, 0.028013])

    # Surface plots for density
    plt.figure(figsize=(8, 6))
    cp_density = plt.contourf(P_grid / 1e6, T_grid, density_grid, cmap='coolwarm')
    plt.title('Mixture Density')
    plt.xlabel('Pressure (MPa)')
    plt.ylabel('Temperature (K)')
    cbar_density = plt.colorbar(cp_density, label='Density (kg/m^3)')
    plt.show()

    # Surface plots for viscosity
    plt.figure(figsize=(8, 6))
    cp_viscosity = plt.contourf(P_grid / 1e6, T_grid, viscosity_grid, cmap='coolwarm')
    plt.title('Mixture viscosity')
    plt.xlabel('Pressure (MPa)')
    plt.ylabel('Temperature (K)')
    cbar_viscosity = plt.colorbar(cp_viscosity, label='Viscosity (Pa·s)')
    plt.show()


