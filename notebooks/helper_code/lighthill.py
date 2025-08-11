import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy import signal


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

                    # For the moment let's skip summing values inside the conduit. 
                    if Y[j] < 0:
                        continue

                    if t_ret >= 0 and t_ret < t_range[-1]:
                        interpolator = RegularGridInterpolator(points, source, method='linear', bounds_error=False, fill_value=0)

                        src_val = interpolator([t_ret, Y[j], np.linalg.norm([X[i], Z[k]])])

                        if not np.isnan(src_val):
                            p_t[idx] += (src_val / (4 * np.pi * r) ) * dv
                            rs.append(r)
                            src_values.append(src_val)


    print(f"Number of contributions: {len(rs)}") 
    print(f"r average: {np.average(rs)}")
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

            r = np.linalg.norm(np.array([x, y, x]) - x_obs)

            t_ret = t_range[t_idx] - r / c0

            if t_ret > 0 and t_ret < t_range[-1]:

                interpolator = RegularGridInterpolator(points, div_T, method='linear', bounds_error=False, fill_value=0)

                src_val_0 = interpolator([0, t_ret, z, np.linalg.norm([x, y])])
                src_val_1 = interpolator([1, t_ret, z, np.linalg.norm([x, y])])

                if not np.isnan(src_val_0) and not np.isnan(src_val_1):
                    p_t += 1/a * (src_val_0*z + src_val_1*np.linalg.norm([x, y])) / (4 * np.pi * r) * ds
    
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