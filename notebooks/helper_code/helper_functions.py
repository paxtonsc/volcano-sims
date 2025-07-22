import os
from slip_imports import readwritedatafiles
import numpy as np

def get_local_solver_from_index_func(folder, file_prefix, base_path='/Users/paxton/git/volcano_sims'):
    solver_from_i = lambda i : readwritedatafiles.read_data_file(f"{base_path}/{folder}/{file_prefix}_{i}.pkl")

    return solver_from_i

def get_quantities_at_conduit_exit(solver_func, iterations=100, R=10, compute_temp=False, conduit_index=-1):
    u_vec = []
    t_vec = []
    p_vec = []
    slip_vec = []
    temp_vec = []

    for i in range(0, iterations, 1):
        solver = solver_func(i)
        momentum = solver.state_coeffs[:,:,solver.physics.get_momentum_slice()]
        rhoSlip = solver.state_coeffs[:, :, solver.physics.get_state_index("rhoSlip")]
        pressure = solver.physics.compute_additional_variable("Pressure", solver.state_coeffs, True)
        rho = np.sum(solver.state_coeffs[:, :, solver.physics.get_mass_slice()],axis=2,keepdims=True)
        
        temp = np.zeros_like(pressure)
        if compute_temp:
            temp = solver.physics.compute_additional_variable("Temperature", solver.state_coeffs, True)

        # Define velocity as momentum divided by density
        u = momentum.ravel() / rho.ravel()
        slip = rhoSlip.ravel() / rho.ravel()

        # Take only the exit velocity
        u_vec.append(u[conduit_index])
        t_vec.append(solver.time)
        p_vec.append(pressure[conduit_index].ravel())
        slip_vec.append(slip[conduit_index].ravel())
        temp_vec.append(temp[conduit_index].ravel())

    return np.array(t_vec), np.array(p_vec), np.array(slip_vec), np.array(u_vec), np.array(temp_vec)


def get_quantities_at_all_space(solver_func, iter, R=10):
    solver = solver_func(iter)
    
    iarhoA, iarhoWv, iarhoM, imom, ie, iarhoWt, iarhoC, iarhoFm, irhoslip = solver.physics.get_state_indices()

    rhoA = solver.state_coeffs[:, :, iarhoA]
    rhoWv = solver.state_coeffs[:, :, iarhoWv]
    rhoM = solver.state_coeffs[:, :, iarhoM]
    momentum = solver.state_coeffs[:, :, imom]
    energy = solver.state_coeffs[:, :, ie]
    rhoWt = solver.state_coeffs[:, :, iarhoWt]
    rhoC = solver.state_coeffs[:, :, iarhoC]
    rhoF = solver.state_coeffs[:, :, iarhoFm]
    rhoSlip = solver.state_coeffs[:, :, irhoslip]

    return rhoA, rhoWv, rhoM, momentum, energy, rhoWt, rhoC, rhoF, rhoSlip
 