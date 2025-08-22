import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy import signal
from .coordinate_transforms import xyz_to_rz, rz_to_xyz


def compute_monopole_dipole_flux(x_obs, t, x_prime, dS, n, c0, t_f, rho_func, v_func, sigma_func, dt=0.05, dx=0.1):
    """
    x_obs - Observation point (2D array of shape (2,))
    t - Time at which to compute the pressure
    x_prime - Source points (2D array of shape (N, 2))
    dS - Surface element area at each source point (1D array of shape (N,))
    n - Normal vector at each source point (2D array of shape (N, 2))
    c0 - Speed of sound
    t_f - Final time
    rho_func - lambda function that outputs rho given [t, r, z]
    v_func - lambda function that outputs v = [v_r, v_z] given [t, r, z]
    sigma_func - lambda function that outputs stress tensor of shape [2, 2] given [t, r, z]
    points - Points for interpolation grid
    dt - Time step for interpolation (default: 0.01)
    dx - Spatial step for interpolation (default: 0.01)m.
    """

    N = len(x_prime)

    M_mon_dt_vec = np.zeros(N)
    D_div_vec = np.zeros(N)

    for k in range(N):
        # Compute the distance from the source to the observation point
        r = np.linalg.norm(x_obs - x_prime[k])

        x_prime_rz = xyz_to_rz(x_prime[k])
        n_rz = xyz_to_rz(n[k])
        
        # Compute the time of arrival at the observation point
        t_ret = t - r / c0
        
        if t_ret < 0 or t_ret + dt > t_f:
            continue

        #print(f"made it past : {t} at k {k}")
        
        # Interpolate the density and velocity at t=t_ret and x_prime[k]
        rho = rho_func(t_ret, x_prime_rz[0], x_prime_rz[1])
        v = v_func(t_ret, x_prime_rz[0], x_prime_rz[1])

        # Interpolate the density and velocity at t=t_ret + dt and x_prime[k]
        rho_plus = rho_func(t_ret + dt, x_prime_rz[0], x_prime_rz[1])
        v_plus = np.nan_to_num(v_func(t_ret + dt, x_prime_rz[0], x_prime_rz[1]))

        # Calculate the time derivative of the monopole moment.
        # Do calculations in spherical coordinates
        rho_v_n = rho * np.dot(v, n_rz)
        rho_v_n_plus = rho_plus * np.dot(v_plus, n_rz)

        # Calculate the monopole moment time derivative
        M_mon_dt_vec[k] = (rho_v_n_plus - rho_v_n) / dt / (4 * r * np.pi) * dS[k]

        # Calculate the divergence of the dipole term.
        for i in range(3):  # Loop over x and y dimensions for 2D divergence
            x_plus = x_prime[k].copy()
            x_plus[i] += dx
            x_plus_rz = xyz_to_rz(x_plus)

            r_plus = np.linalg.norm(x_obs - x_plus)
            t_ret_plus = t - r_plus / c0

            x_minus = x_prime[k].copy()
            x_minus[i] -= dx
            x_minus_rz = xyz_to_rz(x_minus)

            r_minus = np.linalg.norm(x_obs - x_minus)
            t_ret_minus = t - r_minus / c0

            # Interpolate in polar coordinates.
            rho = rho_func(t_ret_plus, x_plus_rz[0], x_plus_rz[1])
            sigma = sigma_func(t_ret_plus, x_plus_rz[0], x_plus_rz[1])
            v_rz = v_func(t_ret_plus, x_plus_rz[0], x_plus_rz[1])

            # Convert velocity from cylindrical -> cartesian
            v = rz_to_xyz(v_rz.flatten(), n[k])

            # Compute force term: n[k] · (ρ v v + σ) / (r_plus)
            force_plus = np.sum(n[k] * (rho * np.outer(v, v) + sigma), axis=1) / r_plus

            rho = rho_func(t_ret_minus, x_minus_rz[0], x_minus_rz[1])
            sigma = sigma_func(t_ret_minus, x_minus_rz[0], x_minus_rz[1])
            v_rz = v_func(t_ret_minus, x_minus_rz[0], x_minus_rz[1])
            
            # Convert velocity from cylindrical -> cartesian
            v = rz_to_xyz(v_rz.flatten(), n[k])

            # Compute force term: n[k] · (ρ v v + σ) / (r_minus)
            force_minus = np.sum(n[k] * (rho * np.outer(v, v) + sigma), axis=1) / r_minus

            # Accumulate divergence for component i
            D_div_vec[k] += (force_plus[i] - force_minus[i]) / (2 * dx) / (4 * np.pi) * dS[k]

    return np.sum(M_mon_dt_vec), np.sum(D_div_vec), M_mon_dt_vec, D_div_vec


def compute_ref_mapping(x_tri):
	''' Compute jacobian
	x_tri: triangle nodes (..., 3, 2)

	Returns
	Matrix collection (..., 2, 2) such that reference coordinates can be
	obtained from

		M @ (x - x_tri[0,:]).

	'''
	_a00 = x_tri[...,1:2,0:1] - x_tri[...,0:1,0:1]
	_a01 = x_tri[...,2:3,0:1] - x_tri[...,0:1,0:1]
	_a10 = x_tri[...,1:2,1:2] - x_tri[...,0:1,1:2]
	_a11 = x_tri[...,2:3,1:2] - x_tri[...,0:1,1:2]
	_M = np.zeros((*x_tri.shape[:-2], 2, 2))
	_M[...,0:1,0:1] = _a11
	_M[...,0:1,1:2] = -_a01
	_M[...,1:2,0:1] = -_a10
	_M[...,1:2,1:2] = _a00
	_M /= (_a00 * _a11 - _a01 * _a10)
	return _M



def calculate_pressure_as_volume_integral(X, Y, Z, file_index_list, x_obs, points, source, t_range, max_r, c0):
    dv = (X[1] - X[0]) * (Y[1] - Y[0]) * (Z[1] - Z[0])

    print(f"DV size is {dv}")
    print(f"X size is {len(X)}, Y size is {len(Y)}, Z size is {len(Z)}, t_range size is {len(t_range)}, file_index_list size is {len(file_index_list)}, x_obs size is {len(x_obs)}, points size is {len(points)}")

    rs = []
    src_values = []
    p_t = np.zeros(len(file_index_list))
    src_values = []

    for idx, t_idx in enumerate(file_index_list):
        for i in range(len(X)):
            for  j in range(len(Y)):
                for k in range(len(Z)):
                    # source position
                    y_src = np.array([X[i], Y[j], Z[k]])

                    # distance from source to observation point
                    r = np.linalg.norm(x_obs - y_src)
                    # change?
                    t_ret = t_range[t_idx] - r / c0

                    # a little hacky, but ignore any contributions with r the max value of Y
                    if np.linalg.norm([X[i], Z[k]]) > max_r:
                        continue

                    if t_ret >= 0 and t_ret < t_range[-1]:
                        interpolator = RegularGridInterpolator(points, source, method='linear', bounds_error=False, fill_value=0)

                        src_val = interpolator([t_ret, Y[j], np.linalg.norm([X[i], Z[k]])])
                        src_values.append(src_val)

                        if not np.isnan(src_val):
                            p_t[idx] += (src_val / (4 * np.pi * r) ) * dv
                            rs.append(r)
                            src_values.append(src_val)


    print(f"Number of contributions: {len(rs)}") 
    print(f"r average: {np.average(rs)}")
    print(f"source value average: {np.average(np.asarray(src_values))}")
    #print(f"Max source value: {max(np.asarray(src_values))}")
    return p_t


def calculate_surface_integral(div_T, a, t_idx, t_range, x_obs, c0, points, N_theta=50, N_phi=50):

    theta = np.linspace(0, np.pi/3, N_theta)
    phi = np.linspace(0, 2*np.pi, N_phi)

    dphi = phi[1] - phi[0]
    dtheta = theta[1] - theta[0]

    p_t = 0

    for i in range(N_theta):
        for j in range(N_phi):
            ds = a**2 * np.sin(phi[j]) * dphi * dtheta
            x = a * np.sin(phi[j]) * np.cos(theta[i])
            y = a * np.sin(phi[j]) * np.sin(theta[i])
            z = a * np.cos(phi[j])

            r = np.linalg.norm(np.array([x, y, z]) - x_obs)

            t_ret = t_range[t_idx] - r / c0

            if t_ret > 0 and t_ret < t_range[-1]:

                interpolator = RegularGridInterpolator(points, div_T, method='linear', bounds_error=False, fill_value=0)

                src_val_0 = interpolator([0, t_ret, z, np.linalg.norm([x, y])])
                src_val_1 = interpolator([1, t_ret, z, np.linalg.norm([x, y])])

                if not np.isnan(src_val_0) and not np.isnan(src_val_1):
                    p_t += 1/a * (src_val_0*z + src_val_1*np.linalg.norm([x, y])) / (4 * np.pi * r) * ds
    
    return p_t


def mario_surface_integral(a, t_idx, t_range, x_obs, c0, rho, v_r, v_z, sigma, points, N_theta=50, N_phi=50):

    theta = np.linspace(0, np.pi/3, N_theta)
    phi = np.linspace(0, 2*np.pi, N_phi)

    dphi = phi[1] - phi[0]
    dtheta = theta[1] - theta[0]

    p_t = 0

    for i in range(N_theta):
        for j in range(N_phi):
            ds = a**2 * np.sin(phi[j]) * dphi * dtheta
            x = a * np.sin(phi[j]) * np.cos(theta[i])
            y = a * np.sin(phi[j]) * np.sin(theta[i])
            z = a * np.cos(phi[j])

            r = np.linalg.norm(np.array([x, y, z]) - x_obs)

            t_ret = t_range[t_idx] - r / c0

            if t_ret > 0 and t_ret < t_range[-1]:

                v_r_interpolator = RegularGridInterpolator(points, v_r, method='linear', bounds_error=False, fill_value=0)
                v_z_interpolator = RegularGridInterpolator(points, v_z, method='linear', bounds_error=False, fill_value=0)
                sigma_interpolator = RegularGridInterpolator(points, sigma, method='linear', bounds_error=False, fill_value=0)
                rho_interpolator = RegularGridInterpolator(points, rho, method='linear', bounds_error=False, fill_value=0)

                v_r = v_r_interpolator([t_ret, z, np.linalg.norm([x, y])])
                v_z = v_z_interpolator([t_ret, z, np.linalg.norm([x, y])])
                sigma = sigma_interpolator([t_ret, z, np.linalg.norm([x, y])])
                rho = rho_interpolator([t_ret, z, np.linalg.norm([x, y])])

                p_t += 1/a * (v_r*z + v_z*np.linalg.norm([x, y])) / (4 * np.pi * r) * ds

    return p_t

def highpass(p, lowcut=1.0):
    fs = 10  # Sampling frequency (Hz)

    # Step 2: Define band-pass filter parameters
    lowcut = 1.0  # Lower cutoff frequency (Hz)

    order = 5  # Filter order

    # Normalize cutoff frequencies to Nyquist frequency (fs/2)
    nyquist = 0.5 * fs
    low = lowcut / nyquist

    # Step 3: Design Butterworth band-pass filter
    b, a = signal.butter(order, low, btype='high')

    return signal.filtfilt(b, a, p)

def lowpass_weighted_average(p, point_distance=1):

    averaged_p = np.zeros_like(p)

    for i in range(len(p)):
        start = max(0, i - point_distance)
        end = min(len(p), i + point_distance + 1)
        averaged_p[i] = np.mean(p[start:end])

    return averaged_p