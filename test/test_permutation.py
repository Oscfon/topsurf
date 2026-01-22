import pytest

def test_str_to_cycles():
    from topsurf.permutation import str_to_cycles

    assert str_to_cycles("(0,1)") == [[0, 1]]
    assert str_to_cycles("(0,1)(3,2)") == [[0, 1], [3, 2]]
    assert str_to_cycles("()(0,1)()(2,3)") == [[0, 1], [2, 3]]
    assert str_to_cycles("(0,1,2)(~0,~1,~2)") == [[0, 1, 2], [-1, -2, -3]]

    with pytest.raises(TypeError):
        str_to_cycles(2)


def test_order():
    from topsurf.permutation import perm_init, perm_order

    p = perm_init("()")
    assert perm_order(p) == 1

    p = perm_init("(1)")
    assert perm_order(p) == 1

    p = perm_init("(1,2)")
    assert perm_order(p) == 2

    p = perm_init("(1,2)(3,4)(6,7,8)")
    assert perm_order(p) == 6
