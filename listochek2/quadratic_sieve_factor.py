import math

import libnum
import numpy.linalg
import numpy as np

primes_b = []
a_splits = []


def modular_square_root(n, p, k):
    if libnum.has_sqrtmod_prime_power(n, p, k):
        return list(libnum.sqrtmod_prime_power(n, p, k))
    return []


# returns list of pairs (p, k) where p is prime and p ** k <= n
def prime_powers(b, n):
    primes = []
    prime = [True] * (b + 1)
    prime[0] = prime[1] = False

    for i in range(2, b):
        if prime[i]:
            primes_b.append(i)
            j = 1
            while i ** j < n:
                primes.append((i, j))
                j += 1
            for j in range(i * i, b + 1, i):
                prime[j] = False

    for i in range(int(n ** 0.5) + 1, b + 1):
        if prime[i]:
            primes.append((i, 1))

    return primes


def quadratic_sieve(b, n):
    factor_base = prime_powers(b, n)

    def f(t):
        return t ** 2 - n

    aa = int(math.sqrt(n)) + 1
    curr = aa
    fs = [f(curr)]
    while True:
        if fs[-1] > n:  # replace 3n on n
            break
        curr += 1
        fs.append(f(curr))

    for pair in factor_base:
        root = modular_square_root(n, pair[0], pair[1])
        step = pair[0] ** pair[1]
        for el in root:
            for i in range((el + step - aa % step) % step, len(fs), step):
                fs[i] //= pair[0]

    smooth_squares = []
    for i in range(len(fs)):
        if fs[i] == 1:
            smooth_squares.append(i + aa)

    return smooth_squares


def split_c(c):
    i = 0
    k = [0] * len(primes_b)
    while c > 1:
        if c % primes_b[i] == 0:
            k[i] += 1
            c //= primes_b[i]
        else:
            i += 1
    a_splits.append(k.copy())
    for j in range(len(k)):
        k[j] %= 2
    return k


def quadratic_sieve_factor(N):
    global a
    a = quadratic_sieve(100, N)
    k = [split_c(el * el % N) for el in a]
    k = numpy.array(k).transpose()
    rows, cols = k.shape
    max_size = max(rows, cols)
    square_k = np.zeros((max_size, max_size), dtype=k.dtype)
    square_k[:rows, :cols] = k
    solve(square_k)
    return result


def solve(k):
    n = k[0].size
    for j in range(n):
        if k[j][j] != 1:
            for i in range(j + 1, n):
                if k[i][j] == 1:
                    k[j] = (k[i] + k[j]) % 2
                    break
        for i in range(j + 1, n):
            if k[i][j] == 1:
                k[i] = (k[j] + k[i]) % 2
    search(k, np.array([0] * n), n - 1)


def search(k, v, i):
    global result
    if result != 1:
        return
    if i == -1:
        check(v)
        return
    s = 0
    for j in range(i + 1, v.size):
        s += v[j] * k[i][j]
    s %= 2
    if k[i][i] == 0:
        if s == 0:
            v[i] = 0
            search(k, v, i - 1)
            v[i] = 1
            search(k, v, i - 1)
        else:
            return
    else:
        if s == 0:
            v[i] = 0
            search(k, v, i - 1)
        else:
            v[i] = 1
            search(k, v, i - 1)


def check(v):
    global a, a_splits, N, result

    prod = 1
    for i in range(len(a)):
        if v[i] == 1:
            prod *= a[i]

    square = 1
    for i in range(len(a)):
        if v[i] == 1:
            for j in range(len(primes_b)):
                if a_splits[i][j] != 0:
                    square *= primes_b[j] ** a_splits[i][j]
    gcd = math.gcd(N, prod - int(math.sqrt(square)))
    if gcd != 1 and gcd != N:
        result = gcd


result = 1
N = 9788111
print(quadratic_sieve_factor(N))
