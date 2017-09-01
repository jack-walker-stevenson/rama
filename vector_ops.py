import operator
from math import sqrt, atan2, degrees


# pythonic vector ops are extensible to higher dimensions but slower here
def add(v0, v1):
    return v0[0] + v1[0], v0[1] + v1[1], v0[2] + v1[2]


def add_pythonic(v0, v1):
    return tuple(map(operator.add, v0, v1))


def subtract(v0, v1):
    return v0[0] - v1[0], v0[1] - v1[1], v0[2] - v1[2]


def subtract_pythonic(v0, v1):
    return tuple(map(operator.sub, v0, v1))


def normalize(v):
    mag = sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    return v[0] / mag, v[1] / mag, v[2] / mag


def normalize_pythonic(v):
    mag = sqrt(sum([x ** 2 for x in v]))
    return tuple([x / mag for x in v])


def dot(v0, v1):
    return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2]


def dot_pythonic(v0, v1):
    return sum([pair[0] * pair[1] for pair in zip(v0, v1)])


def cross(v0, v1):
    return (v0[1] * v1[2] - v0[2] * v1[1], -(v0[0] * v1[2] - v0[2] * v1[0]),
            v0[0] * v1[1] - v0[1] * v1[0])


def torsion(a, b, c, d):
    b1 = subtract(b, a)
    b2 = subtract(c, b)
    b3 = subtract(d, c)
    n1 = normalize(cross(b1, b2))
    n2 = normalize(cross(b2, b3))
    m1 = cross(n1, normalize(b2))
    x = dot(n1, n2)
    y = dot(m1, n2)
    return -1 * degrees(atan2(y, x))
