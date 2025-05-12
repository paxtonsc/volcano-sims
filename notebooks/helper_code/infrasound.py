import numpy as np
from scipy.interpolate import interp1d

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


def find_elem_ID(x, y, trifinders):
  ''' Returns element ID corresponding to given x and y. '''

  trifinders_at_point = [trifinders[0](x, y), trifinders[1](x, y), trifinders[2](x, y)]

  elem_ID = max(trifinders_at_point)
  index = trifinders_at_point.index(elem_ID)

  return index, elem_ID


def get_p_series(x_target, y_target, solver2D_1, solver2D_2, solver2D_3, trifinders, iterations=100):

  # Find the element ID
  index, elem_ID = find_elem_ID(x_target, y_target, trifinders)

  print(f"Element ID for point ({x_target}, {y_target}): {elem_ID}")

  p_relative_arr = []
  p0 = None

  for i in range(0, iterations):
      # Get the solver state at the current time step
      solvers = [solver2D_1(i), solver2D_2(i), solver2D_3(i)]

      U = solvers[index].state_coeffs
      U_target = U[elem_ID:elem_ID+1]

      # Compute pressure using the state vector
      p_target = solvers[index].physics.compute_variable("Pressure", U_target)

      if i == 0:
          p0 = np.average(p_target)

      p_relative_arr.append(np.average(p_target) - p0)
  
  return p_relative_arr