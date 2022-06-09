


import pyMOE as moe

milli = 1e-3
micro = 1e-6
nano = 1e-9
N = 50

def test_GDSMask():

    mask = moe.generate.create_empty_aperture(-500*micro, 500*micro, N, -500*micro, 500*micro, N,)
    # f=50mm, lambda=532nm, R=500µm
    mask = moe.generate.fresnel_phase(mask, 50*milli, 532*nano, radius=500*micro)
    mask.discretize(4)
    moe.plotting.plot_aperture(mask, )
    gdsmask = moe.GDSMask(mask)


def test_GDSMask_create_layout():

    mask = moe.generate.create_empty_aperture(-500*micro, 500*micro, N, -500*micro, 500*micro, N,)
    # f=50mm, lambda=532nm, R=500µm
    mask = moe.generate.fresnel_phase(mask, 50*milli, 532*nano, radius=500*micro)
    mask.discretize(4)
    moe.plotting.plot_aperture(mask, )
    gdsmask = moe.GDSMask(mask)

    gdsmask.create_layout(merge=True)