from hypothesis import given
from hypothesis.strategies import integers, tuples, floats, nothing
from hypothesis.extra.numpy import arrays

import numpy as np
from scipy.spatial.distance import cdist

from pdbear.chirality import set_zero_point

@given(
        values=arrays(
            np.float64, 
            shape=tuples(integers(min_value=1, max_value=10), integers(3,3)),
            elements=floats(min_value=-10, max_value=10),
            fill=nothing(),
            ),
        zero_point=arrays(
            np.float64,
            shape=3,
            elements=floats(min_value=-10, max_value=10),
            fill=nothing(),
            )
        )
def test_distance_remain_same(values, zero_point):
    """Test wether distances between all points remain the same"""

    distances = cdist(values, values)
    new_values = set_zero_point(values, zero_point)
    new_distances = cdist(new_values, new_values)
    assert np.allclose(distances, new_distances)
