from fractions import Fraction
import itertools

from sympy import zeros


def norm_and_trace(coeffs, f):
    n = len(coeffs)

    # Матрица умножения на элемент
    M = [[Fraction(0) for _ in range(n)] for _ in range(n)]

    for i in range(n):
        M[0][i] = Fraction(coeffs[i])

    for i in range(1, n):
        M[i] = polynomial_remainder([Fraction(0)] + M[i - 1], f)
    # for row in M:
    #     print(list(map(str, row)))

    trace = sum(M[i][i] for i in range(n))

    norm = determinant(M)

    return norm, trace


def determinant(matrix):
    n = len(matrix)

    # Базовый случай для матрицы 1x1
    if n == 1:
        return matrix[0][0]

    det = Fraction(0)
    for j in range(n):
        submatrix = [row[:j] + row[j + 1:] for row in matrix[1:]]
        det += (-1) ** j * matrix[0][j] * determinant(submatrix)

    return det


def polynomial_remainder(A, B):
    n = len(A) - 1
    m = len(B) - 1

    if m > n:
        return A

    remainder = A[:]

    while len(remainder) - 1 >= m:
        coeff = remainder[-1] / B[-1]

        temp = [coeff * b for b in B] + [Fraction(0)] * (len(remainder) - len(B))

        remainder = [r - t for r, t in zip(remainder + [Fraction(0)] * (len(temp) - len(remainder)), temp)]

        while remainder and abs(remainder[-1]) < 1e-10:
            remainder.pop()

    remainder_length = n
    remainder += [Fraction(0)] * (remainder_length - len(remainder))

    return remainder


def test():
    # f(x) = x^3 - 2 (неприводим над Z)
    f = [Fraction(-2), Fraction(0), Fraction(0), Fraction(1)]

    # Элемент 1 + θ + 2θ^2
    coeffs = [Fraction(1), Fraction(1), Fraction(2)]

    norm, trace = norm_and_trace(coeffs, f)
    print(f"Норма: {norm}") # 23
    print(f"След: {trace}") # 3

def find_algebraic_integers(f, p):
    for lamb in itertools.product(range(p), repeat=len(f) - 1):
        if lamb.count(0) == len(lamb):
            continue
        norm, trace = norm_and_trace([Fraction(x, p) for x in lamb], f)
        if int(trace) == trace and int(norm) == norm:
            print(lamb, f"Trace:{trace}", f"Norm:{norm}")

# Test 5**(1/3) from the "Computing Integral Bases by John Paul Cook"
f = [Fraction(-5), Fraction(0), Fraction(0), Fraction(1)]
find_algebraic_integers(f, 5)
find_algebraic_integers(f, 3)

# Test 7**(1/3)
f = [Fraction(-7), Fraction(0), Fraction(0), Fraction(1)]
find_algebraic_integers(f, 7)
find_algebraic_integers(f, 3)