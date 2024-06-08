import math

import sympy
import galois
from sympy import Poly

x1, x2 = sympy.symbols("x y")


def phi(m, x, y):
    return x * ksi(m, x, y) ** 2 - ksi(m + 1, x, y) * ksi(m - 1, x, y)


def omega(m, x, y):
    return (ksi(m + 2, x, y) * ksi(m - 1, x, y) ** 2 - ksi(m - 2, x, y) * ksi(m + 1, x, y) ** 2) // (4 * y)


def poly_ksi(n):
    global a, b, q, x1, x2
    ksis = [-1] * (n + 1)

    def f(i):
        if ksis[i] != -1:
            return ksis[i]
        elif i == 0:
            ksis[i] = 0
        elif i == 1:
            ksis[i] = 1
        elif i == 2:
            ksis[i] = Poly(2 * x2, modulus=q)
        elif i == 3:
            ksis[i] = Poly(3 * x1 ** 4 + 6 * int(a) * x1 ** 2 + 12 * int(b) * x1 - int(a) ** 2, modulus=q)
        elif i == 4:
            ksis[i] = Poly(4 * x2 * (
                    x1 ** 6 + 5 * int(a) * x1 ** 4 + 20 * int(b) * x1 ** 3 - 5 * int(a) ** 2 * x1 ** 2 - 4 * int(
                a * b) * x1 - 8 * int(b) ** 2 - int(a) ** 3),
                           modulus=q)
        else:
            m = i // 2
            if i % 2 == 1:
                ksis[i] = Poly(f(m + 2) * f(m) ** 3 - f(m - 1) * f(m + 1) ** 3, modulus=q)
            else:
                ksis[i] = Poly(
                    Poly(f(m) * (f(m + 2) * f(m - 1) ** 2 - f(m - 2) * f(m + 1) ** 2), modulus=q) // Poly(2 * x1,
                                                                                                          modulus=q),
                    modulus=q)
        return ksis[i]

    f(n)
    for el in ksis:
        print(el)
    return f(n)


def ksi(n, x, y):
    global a, b
    ksis = [-1] * (n + 1)

    def f(i):
        if ksis[i] != -1:
            return ksis[i]
        elif i == 0:
            ksis[i] = gf(0)
        elif i == 1:
            ksis[i] = gf(1)
        elif i == 2:
            ksis[i] = 2 * y
        elif i == 3:
            ksis[i] = 3 * x ** 4 + 6 * a * x ** 2 + 12 * b * x - a ** 2
        elif i == 4:
            ksis[i] = 4 * y * (
                    x ** 6 + 5 * a * x ** 4 + 20 * b * x ** 3 - 5 * a ** 2 * x ** 2 - 4 * a * b * x - 8 * b ** 2 - a ** 3)
        else:
            m = i // 2
            if i % 2 == 1:
                ksis[i] = f(m + 2) * f(m) ** 3 - f(m - 1) * f(m + 1) ** 3
            else:
                ksis[i] = f(m) * (f(m + 2) * f(m - 1) ** 2 - f(m - 2) * f(m + 1) ** 2) // (2 * y)
        return ksis[i]

    return f(n)


def mult_n(x, y, n):
    ksi_res = ksi(n, x, y)

    if ksi_res == 0:
        return "O"
    return phi(n, x, y) / ksi_res ** 2, omega(n, x, y) / ksi_res ** 3


def f(n):
    return gf(int((n + q) % q))


def find_points(m):
    global a, b, q, gf
    if m % 2 != 1:
        print("m must be odd")
        return
    result = ["O"]
    ksi = poly_ksi(m)
    for x in Poly(ksi.as_expr().subs({x2 ** 2: x1 ** 3 + int(a) * x1 + int(b)}), modulus=q).ground_roots():
        roots = Poly(x1 ** 2 - (x ** 3 + int(a) * x + int(b)), modulus=q).ground_roots()
        if roots is not None:
            for root in roots:
                result.append([x, root])
    return result


def is_prime(num):
    if num < 2:
        return False
    for i in range(2, int(num ** 0.5) + 1):
        if num % i == 0:
            return False
    return True


def consecutive_primes(n):
    result = []
    i = 2
    while True:
        if is_prime(i):
            result.append(i)
            product = 1
            for prime in result:
                product *= prime
            if product >= n:
                result.pop()
                break
        i += 1
    return result


def solve2(x):
    global q
    return 1 if math.gcd(x ** q - x, x ** 3 + a * x + b) == 1 else 0


def solve_odd(l, x, y):
    global q, gf
    q = l
    ql = q % l
    if ql > l / 2:
        ql -= q
    ksi_l = ksi(l, x, y)
    x_q2 = x ** (q ** 2)
    y_q2 = y ** (q ** 2)
    x_ql, y_ql = mult_n(x, y, ql)
    x_new = (((y_q2 - y_ql) / (x_q2 - x_ql)) ** 2 - x_q2 - x_ql) % ksi_l
    found = False
    for j in range(1, l / 2):
        x_j, y_j = mult_n(x, y, j)
        if x_new % ksi_l == x_j ** q % ksi_l:
            found = True
            break



n = 1000  # Пример значения n
print(consecutive_primes(n))

q = 3623
gf = galois.GF(q)
a = gf(14)
b = gf(19)
print(mult_n(gf(6), gf(730), 2))
print(mult_n(gf(6), gf(730), 4))
print(mult_n(gf(6), gf(730), 8))
print(mult_n(gf(6), gf(730), 16))

a = 14
b = 19
q = 11
gf = galois.GF(q)
a = gf(3)
b = gf(7)
print(ksi(7, gf(5), gf(2)))
print(mult_n(gf(5), gf(2), 5))
print(find_points(5))
