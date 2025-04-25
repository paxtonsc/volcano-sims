import os
from slip_imports import readwritedatafiles

def get_local_solver_from_index_func(folder, file_prefix, base_path='/Users/paxton/git/volcano_sims'):
    solver_from_i = lambda i : readwritedatafiles.read_data_file(f"{base_path}/{folder}/{file_prefix}_{i}.pkl")

    return solver_from_i
