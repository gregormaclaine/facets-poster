x, y, z = sym.symbols('x, y, z')
F = x**3 + x**2 * y + x * y**2 + y**3 + z*(x**2 + x * y + y**2) \
     - z**2*(x + y) - z**3
O = [1, -1, 0]

# Returns whether two points are equal down to scaling
def points_equal(P1, P2):
    return x in sym.solve([p1 - x * p2 for p1, p2 in zip(P1, P2)], x)

# Filters out the two points we know and returns the new point found
def get_new_point(points, given):
    return [p for p in points if not 
any(points_equal(p, g) for g in given)][0]

# Differentiates F and returns the tangent line eq at point P
def tangent_at_point(P):
    a = sym.diff(F, x).subs(x, P[0]).subs(y, P[1]).subs(z, P[2])
    b = sym.diff(F, y).subs(x, P[0]).subs(y, P[1]).subs(z, P[2])
    c = sym.diff(F, z).subs(x, P[0]).subs(y, P[1]).subs(z, P[2])
    return a * x + b * y + c * z

# Takes in two points and performs the algorithm defined in [2]
def add_points(P1, P2):
    L = tangent_at_point(P1) if points_equal(P1, P2) else \
        sym.Matrix([P1, P2, [x, y, z]]).det()
    points_on_L = sym.solve((F, L, x - 1), x, y, z)

    Q = get_new_point(points_on_L, [P1, P2])
    L2 = tangent_at_point(Q) if points_equal(Q, O) else \
        sym.Matrix([Q, O, [x, y, z]]).det()
    points_on_L2 = sym.solve((F, L2, x - 1), x, y, z)

    return integerise(get_new_point(points_on_L2, [Q, O]))


# Take point [1 : 0 : 1] and add it to itself 3 times, creating 8[1 : 0 : 1]
P1 = (1, 0, 1)
_2xP1 = add_points(P1, P1)
_4xP1 = add_points(_2xP1, _2xP1)
print(add_points(_4xP1, _4xP1))