
import pyMOE as moe



def test_create_empty_aperture():
    micro = 1e-6
    mask = moe.generate.create_empty_aperture(-500*micro, 500*micro, 1001, -500*micro, 500*micro, 1001,)
    assert type(mask) is moe.Aperture

def test_circular_aperture():
    micro = 1e-6
    aperture = moe.generate.create_empty_aperture(-500*micro, 500*micro, 1001, -500*micro, 500*micro, 1001,)
    circle = moe.generate.circular_aperture(aperture, radius=250*micro, center=(100*micro, 150*micro))
    