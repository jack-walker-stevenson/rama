import time
import operator
from random import uniform
from math import sqrt


def add_non_pythonic(v0, v1):
    return v0[0] + v1[0], v0[1] + v1[1], v0[2] + v1[2]


def add_map(v0, v1):
    return tuple(map(operator.add, v0, v1))


def normalize_non_pythonic(v):
    mag = sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
    return v[0] / mag, v[1] / mag, v[2] / mag


def normalize(v):
    mag = sqrt(sum([x ** 2 for x in v]))
    return [x / mag for x in v]


def dot_non_pythonic(v0, v1):
    return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2]


def dot(v0, v1):
    return sum([pair[0] * pair[1] for pair in zip(v0, v1)])


def timing(f, n):
    a = (uniform(-1000, 1000), uniform(-1000, 1000), uniform(-1000, 1000))
    b = (uniform(-1000, 1000), uniform(-1000, 1000), uniform(-1000, 1000))
    print(f.__name__)
    r = range(n)
    t1 = time.clock()
    for i in r:
        f(a, b)
        f(a, b)
        f(a, b)
        f(a, b)
        f(a, b)
        f(a, b)
        f(a, b)
        f(a, b)
        f(a, b)
        f(a, b)
    t2 = time.clock()
    print(round(t2-t1, 3))

timing(dot, 100000)
timing(dot_non_pythonic, 100000)
