from libra_toolbox.neutronics.neutron_source import *
from collections.abc import Iterable

def test_get_avg_neutron_rate():

    source = A325_generator_diamond((0, 0, 0), (0, 0, 1))

    assert isinstance(source, Iterable)
    for s in source:
        assert isinstance(s, openmc.IndependentSource)