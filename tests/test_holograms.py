




import pyMOE as moe
import numpy as np

milli = 1e-3
micro = 1e-6
nano = 1e-9
N = 50


def test_algorithm_Gerchberg_Saxton():

    target = np.random.random((128,128))

    # Binary level phase mask
    levels = 4
    iterations = 2
    levels = moe.utils.create_levels(-np.pi, np.pi, levels,)
    phase_mask, errors = moe.holograms.algorithm_Gerchberg_Saxton(target, iterations=iterations, levels=levels)
    far_field = moe.holograms.calculate_phase_farfield(phase_mask)
