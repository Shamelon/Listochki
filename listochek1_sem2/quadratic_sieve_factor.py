import math
import libnum
import numpy.linalg
import numpy as np
from sympy import isprime

FS_LIMIT = 100_000 # limit for quadratic sieve

primes_b = [] # for primes <= b
a_splits = [] # for factoring vectors
smooth_squares = [] # for smooth_squares


def big_l(x: int) -> float:
    # e^sqrt(lnx * lnlnx)
    return math.exp(math.sqrt(math.log(x) * math.log(math.log(x))))


def modular_square_root(n, p, k):
    n %= p ** k
    if libnum.has_sqrtmod_prime_power(n, p, k):
        return list(libnum.sqrtmod_prime_power(n, p, k))
    return []


# returns list of pairs (p, k) where p is prime and p ** k <= n
def prime_powers(b, n):
    global primes_b
    primes = []
    prime = [True] * (b + 1)
    prime[0] = prime[1] = False

    primes_b = []
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


# returns list of B-smooth numbers
def quadratic_sieve(b, n):
    global smooth_squares
    factor_base = prime_powers(b, n)

    def f(t):
        return t ** 2 - n

    aa = int(math.sqrt(n)) + 1
    curr = aa
    fs = [f(curr)]
    while True:
        if fs[-1] > n or len(fs) > FS_LIMIT:
            fs.pop()
            break
        curr += 1
        fs.append(f(curr))

    for p, k in factor_base:
        root = modular_square_root(n, p, k)
        step = p ** k
        for el in root:
            for i in range((el + step - aa % step) % step, len(fs), step):
                assert fs[i] % p == 0, "Bad division"
                fs[i] //= p

    smooth_squares = []
    for i in range(len(fs)):
        if fs[i] == 1:
            smooth_squares.append(i + aa)


# return a vector v[i] = k if v % primes_b[i] ** k == 0
def split_c(c):
    i = 0
    v = [0] * len(primes_b)
    while i < len(primes_b) and c > 1:
        if c % primes_b[i] == 0:
            v[i] += 1
            c //= primes_b[i]
        else:
            i += 1
    a_splits.append(v.copy())
    return v


def quadratic_sieve_factor(n):
    if (sqrt := int(math.sqrt(n))) ** 2 == n:
        if isprime(sqrt):
            return sqrt
    global smooth_squares, a_splits
    # Define B
    b = int(big_l(n) ** (1 / math.sqrt(2))) + 1

    # Find B-smooth numbers
    quadratic_sieve(b, n)

    # Make matrix
    a_splits = []
    k = [
            [el % 2 for el in split_c(el * el % n)] for el in smooth_squares
        ]
    k = numpy.array(k).transpose()

    # Solve the system
    to_upper_triangular(k)
    return generate_solutions(k, find_free_variables(k), n, np.array([0] * k[0].size))


def to_upper_triangular(matrix):
    for i in range(len(matrix)):
        if i >= matrix[0].size:
            return
        for j in range(i + 1, len(matrix)):
            if matrix[j][i] == 1:
                for k in range(i, matrix[j].size):
                    matrix[j][k] ^= matrix[i][k]

# matrix is upper triangular with either 0 or 1 in the cells
def find_free_variables(matrix):
    free_vars = []
    pivot_cols = set()

    for i in range(len(matrix)):
        for j in range(matrix[0].size):
            if matrix[i][j] == 1:
                pivot_cols.add(j)
                break

    for j in range(matrix[0].size):
        if j not in pivot_cols:
            free_vars.append(j)

    return free_vars

# matrix is upper triangular with either 0 or 1 in the cells
# i - the number of the vector in the matrix that is currently being considered
def generate_solutions(matrix, free_vars, n, vector):
    if len(free_vars) == 0:
        return check(vector[::-1], n)

    var = free_vars[0]

    for value in [0, 1]:
        vector[var] = value
        if (res := generate_solutions(matrix, free_vars[1:], n, vector)) != 1:
            return res
    return 1


# v[i] = 1 if smooth_number[i] is included, 0 otherwise
def check(v, n):
    global smooth_squares, a_splits

    prod = 1
    for i in range(len(smooth_squares)):
        if v[i] == 1:
            prod *= smooth_squares[i]

    square = 1
    for i in range(len(smooth_squares)):
        if v[i] == 1:
            for j in range(len(primes_b)):
                if a_splits[i][j] != 0:
                    square *= primes_b[j] ** a_splits[i][j]
    gcd = math.gcd(n, prod - int(math.sqrt(square)))
    if gcd != 1 and gcd != n:
        return gcd
    return 1

def small_tests():
    assert quadratic_sieve_factor(9788111) in [2741, 3571], "the book case"
    assert quadratic_sieve_factor(10004600129) in [100003, 100043], "mid case"
    assert quadratic_sieve_factor(20009200258) in [100003, 100043, 2], "mid case even"

    assert (ans := quadratic_sieve_factor(361)) in [19], "prime square small, given=" + str(ans)
    assert (ans := quadratic_sieve_factor(6859)) in [19], "prime cube small, given=" + str(ans)
    assert (ans := quadratic_sieve_factor(130321)) in [19], "prime 4-th power small, given=" + str(ans)

    assert quadratic_sieve_factor(10008601849) in [100043], "prime square big"
    assert quadratic_sieve_factor(2000000014) in [1000000007, 2], "big difference"

def big_test():
    assert quadratic_sieve_factor(1_004_977_007_034_839) in [1004977, 1000000007], "prime square"

big_test()

