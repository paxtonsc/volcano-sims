import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt



def get_lumped_solution(K=1e9, L_melt=950, tau_p=1e6, tau_r=0, R=10, L_plug=50, g=9.81, rho=2.6e3, mu=5e4, D_c=3, C=600, p0=11.5e6, p_atm=1e5):
    """
    Solves the lumped model differential equations and returns the position and velocity.
    """

    M_plug = R**2 * np.pi * L_plug * rho    # Mass of the plug
    M_eff = M_plug #* (1 + L_melt / (L_plug*2))

    A = np.pi * R**2  # Cross-sectional area

    # Define the shear stress function, currently piecewise linear definition
    def tau(s):
        #return tau_r - (tau_r - tau_p) * np.exp(s / D_c)

        if s < D_c:
            return tau_p - (tau_p - tau_r) * s / D_c
        else:
            return tau_r

    # Define the system of differential equations
    def system(state, t):
        s1, s2 = state  # s1 is position, s2 is velocity
        
        # ds1/dt = s2
        ds1_dt = s2
        
        # ds2/dt = (-A*K*s2)/(M*L_melt) + ((tau_p - tau_s)*2*pi*R*(L_plug - s1)*s1)/M
        ds2_dt =  (-A * K * s1) / (M_eff * L_melt) + \
                (A * (p0 - p_atm) / M_eff) + \
                - (tau(s1) * 2 * np.pi * R * L_plug) / M_eff - \
                4 * np.pi * mu * L_melt * s2 / M_eff \
                - g
        
        return [ds1_dt, ds2_dt]

    # Initial conditions
    s1_0 = 0.0  # Initial position
    s2_0 = 0.0  # Initial velocity
    initial_state = [s1_0, s2_0]

    # Time points
    t_lumped = np.linspace(0, 20, 500)  # Time from 0 to 20 with 500 points

    # Solve the differential equations
    solution = odeint(system, initial_state, t_lumped)

    # Extract solutions
    s1 = solution[:, 0]  # Position
    s2 = solution[:, 1]  # Velocity

    return t_lumped, s1, s2