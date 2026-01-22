import pytest

def test_gcd():
    from topsurf.arith import gcd_int, gcd_list

    assert gcd_int(1, 0) == gcd_int(0, 1) == 1
    assert gcd_int(2, 3) == gcd_int(3, 2) == 1
    assert gcd_int(10, 6) == gcd_int(6, 10) == 2

    assert gcd_list([]) == 0
    assert gcd_list([4, 6, 8, 10]) == 2
    assert gcd_list([2, 2, 3]) == 1


def test_lcm():
    from topsurf.arith import lcm_int, lcm_list

    assert lcm_int(2, 3) == lcm_int(3, 2) == 6
    assert lcm_int(4, 6) == lcm_int(6, 4) == 12

    assert lcm_list([]) == 1
    assert lcm_list([2, 3, 4]) == 12
