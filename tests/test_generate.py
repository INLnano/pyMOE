
import pyMOE as moe
import numpy as np



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

def test_fresnel_zone_plate():
    micro = 1e-6
    milli = 1e-3
    aperture = moe.generate.create_empty_aperture(-5000*micro, 5000*micro, 1001, -5000*micro, 5000*micro, 1001,)
    mask = moe.generate.fresnel_zone_plate_aperture(aperture, 0.5*milli, 532*micro, radius=5000*micro)


def test_aperture_operand():
    micro = 1e-6
    nano = 1e-9
    milli = 1e-3
    aperture1 = moe.generate.create_empty_aperture(-500*micro, 1500*micro, 1001, -500*micro, 500*micro, 1001,)
    aperture1 = moe.generate.fresnel_phase(aperture1, 50*milli, 532*nano, radius=500*micro)

    aperture2 = moe.generate.create_empty_aperture(-500*micro, 1500*micro, 1001, -500*micro, 500*micro, 1001,)
    aperture2 =  moe.generate.arbitrary_aperture_function(aperture2, moe.sag.spiral, L=8)

    operand = np.multiply

    aperture3 = moe.generate.aperture_operation(aperture1, aperture2, operand)

    assert np.all(operand(aperture1.aperture, aperture2.aperture) == aperture3.aperture)
    
def test_aperture_add():
    micro = 1e-6
    nano = 1e-9
    milli = 1e-3
    aperture1 = moe.generate.create_empty_aperture(-500*micro, 1500*micro, 1001, -500*micro, 500*micro, 1001,)
    aperture1 = moe.generate.fresnel_phase(aperture1, 50*milli, 532*nano, radius=500*micro)

    aperture2 = moe.generate.create_empty_aperture(-500*micro, 1500*micro, 1001, -500*micro, 500*micro, 1001,)
    aperture2 =  moe.generate.arbitrary_aperture_function(aperture2, moe.sag.spiral, L=8)

    aperture3 = moe.generate.aperture_add(aperture1, aperture2)

    assert np.all(np.add(aperture1.aperture, aperture2.aperture) == aperture3.aperture)
        
def test_aperture_subtract():
    micro = 1e-6
    nano = 1e-9
    milli = 1e-3
    aperture1 = moe.generate.create_empty_aperture(-500*micro, 1500*micro, 1001, -500*micro, 500*micro, 1001,)
    aperture1 = moe.generate.fresnel_phase(aperture1, 50*milli, 532*nano, radius=500*micro)

    aperture2 = moe.generate.create_empty_aperture(-500*micro, 1500*micro, 1001, -500*micro, 500*micro, 1001,)
    aperture2 =  moe.generate.arbitrary_aperture_function(aperture2, moe.sag.spiral, L=8)

    aperture3 = moe.generate.aperture_subtract(aperture1, aperture2)

    assert np.all(np.subtract(aperture1.aperture, aperture2.aperture) == aperture3.aperture)
    
def test_aperture_multiply():
    micro = 1e-6
    nano = 1e-9
    milli = 1e-3
    aperture1 = moe.generate.create_empty_aperture(-500*micro, 1500*micro, 1001, -500*micro, 500*micro, 1001,)
    aperture1 = moe.generate.fresnel_phase(aperture1, 50*milli, 532*nano, radius=500*micro)

    aperture2 = moe.generate.create_empty_aperture(-500*micro, 1500*micro, 1001, -500*micro, 500*micro, 1001,)
    aperture2 =  moe.generate.arbitrary_aperture_function(aperture2, moe.sag.spiral, L=8)

    aperture3 = moe.generate.aperture_multiply(aperture1, aperture2)

    assert np.all(np.multiply(aperture1.aperture, aperture2.aperture) == aperture3.aperture)
    