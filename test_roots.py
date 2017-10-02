import roots

def test_quadroots_result():
    assert roots.quad_roots(1.0, 1.0, -12.0) == ((3+0j), (-4+0j))

def test_quadroots_types():
    try:
        roots.quad_roots("", "green", "hi")
    except TypeError as err:
        assert(type(err) == TypeError)

def test_quadroots_zerocoeff():
    try:
        roots.quad_roots(a=0.0)
    except ValueError as err:
        assert(type(err) == ValueError)
