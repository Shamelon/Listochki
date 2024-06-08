import secrets
import time
from math import sqrt
from math import log2
from math import gcd

import galois
import numpy as np


def is_prime(n):
    if n <= 1:
        return False
    if n <= 3:
        return True

    if n % 2 == 0 or n % 3 == 0:
        return False

    for i in range(5, int(sqrt(n)), 6):
        if n % i == 0 or n % (i + 2) == 0:
            return False
    return True


def find_prime_factors(s, n):
    while n % 2 == 0:
        s.add(2)
        n = n // 2

    for i in range(3, int(sqrt(n)), 2):
        while n % i == 0:
            s.add(i)
            n = n // i

    if n > 2:
        s.add(n)


def is_primitive(p, q):
    if not is_prime(p):
        print("p is not prime")
        return False

    s = set()

    phi = p - 1
    find_prime_factors(s, phi)

    gf = galois.GF(p)
    q = gf(q)
    for el in s:
        if q ** (phi // el) == 1:
            return False
    return True


def brute_force(g, h, p, check_primitive=False):
    if not is_prime(p):
        print("p is not prime")
        return -1

    if check_primitive and not is_primitive(p, g):
        print("q is not primitive")
        return -1

    if h == 1:
        return p - 1

    res = 1
    for i in range(1, p):
        res = (res * g) % p
        if res == h:
            return i


def task4(a, b, p):
    gf = galois.GF(p)
    a = gf(a)
    b = gf(b)

    q = [0] * int(log2(p - 1))
    if np.power(b, (p - 1) // 2) == gf(1):
        q[0] = 0
    else:
        q[0] = 1
    z = b
    x = q[0]
    pow2 = 1
    for i in range(1, len(q)):
        z *= np.power(a, -q[i - 1] * pow2)
        pow2 *= 2
        m = (p - 1) // np.power(2, i + 1)
        if np.power(z, m) == gf(1):
            q[i] = 0
        else:
            q[i] = 1
        x += q[i] * pow2
    return x


def babystep_giantstep(g, h, p):
    gf = galois.GF(p)
    n = int(sqrt(calculate_order_in_finite_field(g, p))) + 1
    list1 = [gf(1)] * (n + 1)
    for i in range(1, n + 1):
        list1[i] = list1[i - 1] * g

    list2 = [gf(h)] * (n + 1)
    for i in range(1, n + 1):
        list2[i] = list2[i - 1] // list1[n]

    list1 = list(map(int, list1))
    list2 = list(map(int, list2))
    i, j = find_intersection(list1, list2)
    return i + j * n


def calculate_order_in_finite_field(a, p):
    if not (isinstance(a, int) and isinstance(p, int) and p > 1):
        return "Invalid input. 'a' and 'p' should be integers with p > 1."

    order = 1
    result = a % p

    while result != 1:
        result = (result * a) % p
        order += 1

        if order > p:
            return "The order is greater than the prime number p."

    return order


def binary_search(arr, target):
    left = 0
    right = len(arr) - 1

    while left <= right:
        mid = (left + right) // 2
        if arr[mid] == target:
            return True, mid
        elif arr[mid] < target:
            left = mid + 1
        else:
            right = mid - 1

    return False, None


def find_intersection(list1, list2):
    set2 = set(list2)

    for i, num in enumerate(list1):
        if set2.__contains__(num):
            return i, list2.index(num)

    return -1, -1


def generate_unique_random_numbers(n, big_n):
    random_numbers = set()
    while len(random_numbers) < n:
        new_random_number = secrets.randbelow(big_n)
        random_numbers.add(new_random_number)

    return list(random_numbers)


def probabilistic_collision(g, h, p):
    gf = galois.GF(p)
    big_n = calculate_order_in_finite_field(g, p)
    n = 3 * int(sqrt(big_n))
    g = gf(g)
    h = gf(h)
    powers1 = generate_unique_random_numbers(n, max(n, big_n))
    powers2 = generate_unique_random_numbers(n, max(n, big_n))
    list1 = [int(g ** x) for x in powers1]
    list2 = [int(h * g ** x) for x in powers2]
    i, j = find_intersection(list1, list2)
    if i < 0:
        print("No intersection")
        return -1
    return (powers1[i] - powers2[j] + big_n) % big_n


def ext_euc(a, b):
    if b == 0:
        return 1, 0
    x2 = 1
    x1 = 0
    y2 = 0
    y1 = 1
    while b > 0:
        q = a // b
        r = a - q * b
        x = x2 - q * x1
        y = y2 - q * y1
        a = b
        b = r
        x2 = x1
        x1 = x
        y2 = y1
        y1 = y
    x = x2
    y = y2
    return x, y


def mod_inverse(g, n):
    g = g % n
    if gcd(g, n) == 1:
        return None
    else:
        x, _ = ext_euc(g, n)
        return x % n


def solve_mod_equation(g, h, n):
    inv_g = mod_inverse(g, n)
    if inv_g is None:
        return "Обратный элемент для g по модулю n не существует"
    else:
        x = inv_g % n
        return x


def next(alpha, beta, x, p, g, h):
    if 1 >= x >= 2 * p / 3:
        return alpha, (beta + 1) % (p - 1), (h * x) % p
    elif x >= p / 3:
        return (2 * alpha) % (p - 1), (2 * beta) % (p - 1), pow(x, 2, p)
    elif x >= 0:
        return (alpha + 1) % (p - 1), beta, (g * x) % p


def pollard(g, h, p):
    x, y, alpha, beta, gamma, delta = 1, 1, 0, 0, 0, 0
    while True:
        alpha, bet, x = next(alpha, beta, x, p, g, h)
        for j in range(2):
            gamma, delta, y = next(gamma, delta, y, p, g, h)
        if x == y:
            break
    u = (alpha - gamma) % (p - 1)
    v = (delta - beta) % (p - 1)
    s, y1 = ext_euc(v, p - 1)
    d = gcd(v, p - 1)
    if d == 1:
        return pow(v, -1, p - 1) * u
    w = (s * u) % (p - 1)
    for k in range(d):
        x = w // d + k * (p - 1) // d
        if pow(g, x, p) == h:
            return x


# Tests
def task4_test():
    # x  = 7
    print(task4(3, 11, 17))


def pollard_test():
    # x = 37869
    print(pollard(19, 24717, 48611))


def test(g, h, p):
    gf = galois.GF(p)

    start = time.time()
    x = brute_force(g, h, p)
    print("brute-force x =", x)
    print("g**x =", int(gf(g) ** x), ", h =", h)
    print("time:", time.time() - start)

    start = time.time()
    x = babystep_giantstep(g, h, p)
    print("babystep-giantstep x =", x)
    print("g**x =", int(gf(g) ** x), ", h =", h)
    print("time:", time.time() - start)

    start = time.time()
    x = probabilistic_collision(g, h, p)
    print("probabilistic collision x =", x)
    print("g**x =", int(gf(g) ** x), ", h =", h)
    print("time:", time.time() - start)

    start = time.time()
    x = pollard(g, h, p)
    print("pollard x =", x)
    print("g**x =", int(gf(g) ** x), ", h =", h)
    print("time:", time.time() - start)


def medium_test():
    # book page 83
    # x = 1159
    g = 9704
    h = 13896
    p = 17389
    test(g, h, p)


def big_test():
    g = 123_456_789
    p = 1_000_000_007
    h = pow(g, 987_654_321, p)
    test(g, h, p)


if __name__ == '__main__':
    big_test()
