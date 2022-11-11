import sys
# sys.path.insert(0,'..')


import pyMOE as moe
import numpy as np

print("moe")
np.random.seed(123)


def test_aperture_create_aperture():
    x = np.linspace(-500,500,101)
    y = np.linspace(-500,500,101)
    mask = moe.Aperture(x,y)

    assert mask.aperture.size >0


def test_aperture_set_phase():
    x = np.linspace(-500,500,101)
    y = np.linspace(-500,500,101)
    mask = moe.Aperture(x,y)
    phase = np.random.random(mask.shape)
    phase = np.ones(mask.shape)
    mask.phase = phase

    phase = mask.phase

    assert np.all(mask.phase == phase)



def test_aperture_create_aperture_from_array():
    phase_mask = np.random.random((128,128))
    pixel_size = 1e-6
    mask = moe.generate.create_aperture_from_array(phase_mask, pixel_size=pixel_size, center=True, )
    


