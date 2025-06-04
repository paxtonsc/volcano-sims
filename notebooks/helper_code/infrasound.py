import numpy as np
from scipy.interpolate import interp1d
import matplotlib.tri as tri

T_0 = 288.15 # [K]
M = 0.02896968 # [kg/mol]
R_0 = 8.314462618 #[J/(mol·K)]
p_0 = 101325 # [Pa]
g = 9.81 # [m/s^2]

SUMMIT_TUNGURAHUA = 5023 # [m]

L = 0.0065    # Temperature lapse rate (K/m)
gamma = 1.4   # Specific heat ratio for air

def pressure():
    return p_0 * np.exp(-M * g * (SUMMIT_TUNGURAHUA) / (R_0 * T_0)) # [Pa]


# (kg/m^3)
def density():
    T = T_0 - L * (SUMMIT_TUNGURAHUA)
    P = pressure()

    return P / (R_0 * T)

def sound_speed():
    return 320 #[m/s] estimated speed of sound in air at 5000m

def point_in_volcano(x, y):
    # Define the volcano shape
    # This is a placeholder function; replace with actual volcano shape logic
    if y < -1/2 * x and y < 1/2 * x:
        return True
    else:
        return False

def point_inside_vertical_conduit(x, y):
    # Define the volcano shape
    # This is a placeholder function; replace with actual volcano shape logic
    if y < 1000 and y > -40 and x > 0 and x < 10:
        return True
    else:
        return False

def Q_dot_func(t, t_vec, Q_dot_vec):
    if t < t_vec[0] or t > t_vec[-1]:
        return 0
    
    # Create interpolator for each grid point
    interpolator = interp1d(t_vec, Q_dot_vec, kind='linear', fill_value="extrapolate")
    
    # Evaluate at time t
    interpolated_Q = interpolator(t)
    
    # Linear interpolation
    return interpolated_Q


def relative_pressure(t, x, y, t_vec, Q_dot_vec):
    r = np.sqrt(x**2 + y**2)
    rho = density()

    # 2/3 comes from the fact that 1/3 of the outward area is solid volcano
    return rho * Q_dot_func(t - r / sound_speed(), t_vec, Q_dot_vec) / ((2/3) * np.pi * 4 * r)

def relative_pressure_vertical_conduit(t, x, y, t_vec, Q_dot_vec):
    r = np.sqrt(x**2 + (y+40)**2)
    rho = density()

    # 2/3 comes from the fact that 1/3 of the outward area is solid volcano
    return rho * Q_dot_func(t - r / sound_speed(), t_vec, Q_dot_vec)


def find_elem_ID(x, y, trifinder):
  ''' Returns element ID corresponding to given x and y. '''

  return trifinder(x, y)


def get_p_series(x_target, y_target, solver2D_1, trifinder, iterations=100, p0=None, d_iter=2):


    # Find the element ID
    elem_ID = find_elem_ID(x_target, y_target, trifinder)

    print(f"Element ID for point ({x_target}, {y_target}): {elem_ID}")

    p_relative_arr = []

    for i in range(0, iterations, d_iter):
        # Get the solver state at the current time step
        solver = solver2D_1(i)

        U = solver.state_coeffs

        x_node_elem = solver.mesh.node_coords[solver.mesh.elem_to_node_IDs[elem_ID,:], :]

        try:
            # Evaluate state_coeff at the exact point (x_target, y_target)
            U_target = np.array([tri.CubicTriInterpolator(
                tri.Triangulation(x_node_elem[:,0],x_node_elem[:,1]), # Create local triangulation using x, y of relevant triangle
                U[elem_ID, :, i])(x_target,y_target) for i in range(U.shape[-1])])
            
            # Pad U_target to the right shape for physics.compute_variable
            U_target = U_target[np.newaxis, np.newaxis, :]
        except:
            U_target = U[elem_ID:elem_ID+1]


        # Compute pressure using the state vector
        p_target = np.average(solver.physics.compute_variable("Pressure", U_target))
        if i == 0 and p0 is None:
            p0 = p_target

        p_relative_arr.append(p_target - p0)

    return p_relative_arr