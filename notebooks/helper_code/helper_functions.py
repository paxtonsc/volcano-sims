import os
from slip_imports import readwritedatafiles
import numpy as np

def get_local_solver_from_index_func(folder, file_prefix, base_path='/Users/paxton/git/volcano_sims'):
    solver_from_i = lambda i : readwritedatafiles.read_data_file(f"{base_path}/{folder}/{file_prefix}_{i}.pkl")

    return solver_from_i

def get_quantities_at_conduit_exit(solver_func, iterations=100, R=10, computer_temp=False):
    u_vec = []
    t_vec = []
    p_vec = []
    temp_vec = []

    for i in range(0, iterations, 1):
        solver = solver_func(i)
        momentum = solver.state_coeffs[:,:,solver.physics.get_momentum_slice()]
        pressure = solver.physics.compute_additional_variable("Pressure", solver.state_coeffs, True)
        rho = np.sum(solver.state_coeffs[:, :, solver.physics.get_mass_slice()],axis=2,keepdims=True)

        temp = np.zeros_like(pressure)
        if computer_temp:
            temp = solver.physics.compute_additional_variable("Temperature", solver.state_coeffs, True)

        # Define velocity as momentum divided by density
        u = momentum.ravel() / rho.ravel()

        # Take only the exit velocity
        u_vec.append(u[-1])
        t_vec.append(solver.time)
        p_vec.append(pressure[-1].ravel())
        temp_vec.append(temp[-1].ravel())
    
    return np.array(t_vec), np.array(p_vec), np.array(u_vec), np.array(temp_vec)
