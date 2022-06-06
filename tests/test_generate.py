
import pyMOE as moe



def test_create_empty_aperture():
    micro = 1e-6
    mask = moe.generate.create_empty_aperture(-500*micro, 500*micro, 1001, -500*micro, 500*micro, 1001,)
    assert type(mask) is moe.Aperture

def test_circular_aperture():
    micro = 1e-6
    aperture = moe.generate.create_empty_aperture(-500*micro, 500*micro, 1001, -500*micro, 500*micro, 1001,)
    mask = moe.generate.circular_aperture(aperture, radius=250*micro, center=(100*micro, 150*micro))
    
def test_rectangular_aperture():
    micro = 1e-6
    aperture = moe.generate.create_empty_aperture(-500*micro, 500*micro, 1001, -500*micro, 500*micro, 1001,)
    mask = moe.generate.rectangular_aperture(aperture, width=250*micro, height=100*micro, center=(100*micro, 150*micro))
    mask = moe.generate.rectangular_aperture(aperture, width=250*micro, height=100*micro, corner=(-100*micro, 150*micro))
    mask = moe.generate.rectangular_aperture(aperture, width=250*micro, height=100*micro)
        
def test_fresnel_phase():
    micro = 1e-6
    milli = 1e-3
    aperture = moe.generate.create_empty_aperture(-5000*micro, 5000*micro, 1001, -5000*micro, 5000*micro, 1001,)
    mask = moe.generate.fresnel_phase(aperture, 2*milli, 532*micro, radius=5000*micro)
