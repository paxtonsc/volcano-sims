import os
from slip_imports import readwritedatafiles

def get_local_solver_from_index_func(folder, file_prefix, base_path='/Users/paxton/git/volcano_sims'):
    folder = "mario_file_review/test_plug_boundary"
    file_prefix = "tungurahua_rad_5_v15_conduit_"
    solver_from_i = lambda i : readwritedatafiles.read_data_file(f"{base_path}/{folder}/{file_prefix}{i}.pkl")

    return solver_from_i
