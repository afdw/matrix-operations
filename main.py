from fractions import Fraction
from math import gcd, prod
from collections import Counter
from functools import reduce

try:
    from math import lcm
except:
    def lcm(a, b):
        return abs(a * b) // gcd(a, b)

try:
    gcd(1, 2, 3)
except:
    old_gcd = gcd
    gcd = lambda *x: reduce(old_gcd, x, 0)

try:
    lcm(1, 2, 3)
except:
    old_lcm = lcm
    lcm = lambda *x: reduce(old_lcm, x, 1)

def sign(number):
    if number > 0:
        return 1
    elif number < 0:
        return -1
    else:
        return 0

def number_to_str(number):
    if type(number) is Fraction:
        return f"{number.numerator}/{number.denominator}" if number.denominator != 1 else str(number.numerator)
    else:
        return str(number)

inf = float("inf")

class Polynomial:
    def __init__(self, data = None):
        if data is not None:
            self.data = data
        while self.data and self.data[-1] == 0:
            self.data.pop()
        self.n = len(self.data) - 1 if self.data else -inf
        self.rational = all(type(self.data[i]) in {int, Fraction} for i in Polynomial.range(self.n))
        self.data = [Fraction(self.data[i]) if self.rational else self.data[i] for i in Polynomial.range(self.n)]

    @staticmethod
    def range(n):
        return range(n + 1 if n != -inf else 0)

    def __repr__(self):
        return " + ".join([f"{number_to_str(self[i])} * x^{i}" for i in Polynomial.range(self.n)]) if self.n != 0 else "0"

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        self.data[key] = value
        self.__init__()

    def __len__(self):
        return self.n

    def __iter__(self):
        return (self[i] for i in range(self.n + 1))

    def __eq__(self, other):
        return self.data == other.data

    def __add__(self, other):
        return Polynomial([(self[i] if i <= self.n else 0) + (other[i] if i <= other.n else 0) for i in Polynomial.range(max(self.n, other.n))])

    def __sub__(self, other):
        return Polynomial([(self[i] if i <= self.n else 0) - (other[i] if i <= other.n else 0) for i in Polynomial.range(max(self.n, other.n))])

    def __mul__(self, other):
        if type(other) is not Polynomial:
            return Polynomial([self[i] * other for i in Polynomial.range(self.n)])
        else:
            return Polynomial([sum((self[j] if j <= self.n else 0) * (other[i - j] if i - j <= other.n else 0) for j in Polynomial.range(i)) for i in Polynomial.range(self.n + other.n)])

    def __truediv__(self, other):
        return Polynomial([self[i] / other for i in Polynomial.range(self.n)])

    def __pow__(self, power):
        assert self.n == self.m
        assert type(power) is int
        assert power >= 0
        result = Polynomial([1])
        for _ in range(power):
            result *= self
        return result

    def __call__(self, x):
        """
        >>> Polynomial([7, -9, 2])(Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))
        (58  54  57
         96  124 138
         141 180 226)
        """
        return sum((x**i * self[i] for i in Polynomial.range(self.n)), x - x)

    def monic(self):
        """
        >>> Polynomial([3, 4, 1]) * Polynomial([7, -9, 2])
        21 * x^0 + 1 * x^1 + -23 * x^2 + -1 * x^3 + 2 * x^4
        """
        if self.n != -inf:
            return self / self.data[-1]
        else:
            return self

class Matrix:
    def __init__(self, data = None):
        if data is not None:
            self.data = data
        self.n = len(self.data)
        self.m = 0 if len(self.data) == 0 else len(self.data[0])
        assert all(len(self.data[i]) == self.m for i in range(self.n))
        if self.m == 0:
            self.n = 0
        self.rational = all(type(self.data[i][j]) in {int, Fraction} for i in range(self.n) for j in range(self.m))
        self.data = [[Fraction(self.data[i][j]) if self.rational else self.data[i][j] for j in range(self.m)] for i in range(self.n)]

    @staticmethod
    def identity(n):
        return Matrix([[1 if i == j else 0 for j in range(n)] for i in range(n)])

    def __repr__(self):
        return "\n".join([("(" if i == 0 else " ") + " ".join([(lambda s: s.ljust(max(len(number_to_str(self.data[k][j])) for k in range(self.n))) if i == self.n - 1 or j != self.m - 1 else s)(number_to_str(self.data[i][j])) for j in range(self.m)]) + (")" if i == self.n - 1 else "") for i in range(self.n)])

    def __getitem__(self, key):
        if type(key) is tuple:
            a, b = key
            if type(a) is not slice and type(b) is not slice:
                return self.data[a][b]
            else:
                if type(a) is not slice:
                    a = slice(a, a + 1)
                if type(b) is not slice:
                    b = slice(b, b + 1)
                return Matrix([[x for x in l[b]] for l in self.data[a]])
        else:
            return self.data[key]

    def __setitem__(self, key, value):
        if type(key) is tuple:
            a, b = key
            if type(a) is not slice and type(b) is not slice:
                self.data[a][b] = value
            else:
                if type(a) is not slice:
                    a = slice(a, a + 1)
                if type(b) is not slice:
                    b = slice(b, b + 1)
                for i_1, i_2 in enumerate(range(self.n)[a]):
                    for j_1, j_2 in enumerate(range(self.m)[b]):
                        self.data[i_2][j_2] = value[i_1][j_1]
        else:
            self.data[key] = value
        self.__init__()

    def __len__(self):
        return self.n

    def __iter__(self):
        return (self[i, j] for i in range(self.n) for j in range(self.m))

    def __eq__(self, other):
        return self.data == other.data

    def __add__(self, other):
        assert other.n == self.n and other.m == self.m
        return Matrix([[self[i, j] + other[i, j] for j in range(self.m)] for i in range(self.n)])

    def __sub__(self, other):
        assert other.n == self.n and other.m == self.m
        return Matrix([[self[i, j] - other[i, j] for j in range(self.m)] for i in range(self.n)])

    def __mul__(self, other):
        """
        >>> Matrix.identity(5)[0:3, 0:4] * Matrix.identity(5)[0:4, 0:5]
        (1 0 0 0 0
         0 1 0 0 0
         0 0 1 0 0)
        """
        if type(other) is not Matrix:
            return Matrix([[self[i, j] * other for j in range(self.m)] for i in range(self.n)])
        else:
            assert self.m == other.n or self.n == 0 or other.m == 0
            return Matrix([[sum(self[i, k] * other[k, j] for k in range(self.m)) for j in range(other.m)] for i in range(self.n)])

    def __truediv__(self, other):
        return Matrix([[self[i, j] / other for j in range(self.m)] for i in range(self.n)])

    def __pow__(self, power):
        assert self.n == self.m
        assert type(power) is int
        if power < 0:
            self = self.inverse()
            assert self is not None
            power = -power
        result = Matrix.identity(self.n)
        for _ in range(power):
            result *= self
        return result

    def transpose(self):
        return Matrix([[self[i, j] for i in range(self.n)] for j in range(self.m)])

    def append_right(self, other):
        assert self.n == other.n or self.m == 0 or other.m == 0
        return Matrix([[self[i, j] if j < self.m else other[i, j - self.m] for j in range(self.m + other.m)] for i in range(max(self.n, other.n))])

    def append_bottom(self, other):
        assert self.m == other.m or self.n == 0 or other.n == 0
        return Matrix([[self[i, j] if i < self.n else other[i - self.n, j] for j in range(max(self.m, other.m))] for i in range(self.n + other.n)])

    def append_bottom_right(self, other):
        return Matrix([[self[i, j] if i < self.n and j < self.m else other[i - self.n, j - self.m] if i >= self.n and j >= self.m else 0 for j in range(self.m + other.m)] for i in range(self.n + other.n)])

    def denominators_lcm(self):
        return lcm(*(x.denominator for x in self)) if self.rational else 1

    def row_reduce(self, *, bound=None, integerize=False):
        """
        >>> Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]]).row_reduce()
        (1 0 0
         0 1 0
         0 0 1)
        >>> Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).row_reduce()
        (1 0 -1
         0 1 2
         0 0 0 )
        >>> Matrix([[1, 2, 3, 1, 0, 0], [4, 5, 6, 0, 1, 0], [7, 8, 9, 0, 0, 1]]).row_reduce()
        (1 0 -1 0 -8/3 5/3
         0 1 2  0 7/3  -4/3
         0 0 0  1 -2   1   )
        >>> Matrix([[1, 2, 3, 1, 0, 0], [4, 5, 6, 0, 1, 0], [7, 8, 9, 0, 0, 1]]).row_reduce(integerize=True)
        (3 0 -3 0 -8 5
         0 3 6  0 7  -4
         0 0 0  1 -2 1 )
        """
        if bound is None:
            bound = self.m
        self = Matrix(self)
        l = 0
        for j in range(bound):
            for i in range(l, self.n):
                if self[i, j] != 0:
                    self[l, :], self[i, :] = self[i, :], self[l, :]
                    self[l, :] /= self[l, j]
                    for k in range(self.n):
                        if k != l:
                            self[k, :] -= self[l, :] * self[k, j]
                    l += 1
                    break
        if integerize:
            self *= self.denominators_lcm()
            for i in range(self.n):
                self[i, :] /= max(gcd(*(x.numerator for x in self[i, :])), 1)
        return self

    def column_reduce(self, *, bound=None, integerize=False):
        return self.transpose().row_reduce(bound=bound, integerize=integerize).transpose()

    def remove_zero_rows(self):
        return Matrix([[self[i, j] for j in range(self.m)] for i in range(self.n) if any(self[i, j] != 0 for j in range(self.m))])

    def remove_zero_columns(self):
        return self.transpose().remove_zero_rows().transpose()

    def range(self):
        """
        >>> Matrix([[1, 0, -3, 0, 2, -8], [0, 1, 5, 0, -1, 4], [0, 0, 0, 1, 7, -9], [0, 0, 0, 0, 0, 0]]).range()
        (1 0 0
         0 1 0
         0 0 1
         0 0 0)
        """
        return self.column_reduce(integerize=True).remove_zero_columns()

    def rang(self):
        return self.range().m

    def null(self):
        """
        >>> Matrix([[1, 0, -3, 0, 2, -8], [0, 1, 5, 0, -1, 4], [0, 0, 0, 1, 7, -9], [0, 0, 0, 0, 0, 0]]).null()
        (3  -2 8
         -5 1  -4
         1  0  0
         0  -7 9
         0  1  0
         0  0  1 )
        """
        reduced = self.append_bottom(Matrix.identity(self.m)).column_reduce(bound=self.n, integerize=True)
        return Matrix([[reduced[i, j] for j in range(self.m) if all(reduced[i, j] == 0 for i in range(self.n))] for i in range(self.n, self.n + self.m)])

    def inverse(self):
        """
        >>> Matrix([[1, 2, 3], [4, 2, 2], [5, 1, 7]]).inverse()
        (-2/7 11/42 1/21
         3/7  4/21  -5/21
         1/7  -3/14 1/7  )
        >>> Matrix([[1, 2], [1, 2]]).inverse()
        """
        assert self.n == self.m
        reduced = self.append_right(Matrix.identity(self.n)).row_reduce()
        return reduced[:, self.n:] if reduced[:, :self.n] == Matrix.identity(self.n) else None

    def image(self, space):
        """
        >>> Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).image(Matrix([[1, 0], [-2, -3], [1, 0]]))
        (2
         5
         8)
        """
        assert space.n == self.m or space.m == 0
        return (self * space).range()

    def inverse_image(self, space):
        """
        >>> Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).inverse_image(Matrix([[6], [15], [24]]))
        (1  0
         -2 -3
         1  0 )
        """
        assert space.n == self.n or space.m == 0
        return self.append_right(space).null()[:self.m, :]

    def find_inverse(self, vector):
        """
        >>> Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).find_inverse(Matrix([[6], [15], [24]]))
        (0
         3
         0)
        >>> Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).find_inverse(Matrix([[1], [0], [0]]))
        """
        assert vector.n == self.n
        null = self.append_right(vector).null()
        return next((null[:self.m, j] / null[self.m, j] * (-1) for j in range(null.m) if null[self.m, j] != 0), None)

    def change_basis(self, in_basis, out_basis = None):
        """
        >>> Matrix([[1, 2, 0], [0, 1, 0], [0, 0, 3]]).change_basis(Matrix([[3, 2, 5], [1, 4, 2], [5, 5, 1]]))
        (247/75 208/75 56/75
         -8/75  13/75  -34/75
         -14/15 4/15   23/15 )
        """
        if out_basis is None:
            out_basis = in_basis
        assert in_basis.inverse() is not None
        assert out_basis.inverse() is not None
        assert in_basis.n == self.m
        assert out_basis.n == self.n
        return out_basis.inverse() * self * in_basis

    def restrict(self, space):
        """
        >>> Matrix([[1, 2, 0], [0, 1, 0], [0, 0, 3]]).restrict(Matrix([[1, 0], [0, 1], [0, 0]]))
        (1 2
         0 1)
        >>> Matrix([[1, 2, 0], [0, 1, 0], [0, 0, 3]]).restrict(Matrix([[1, 0], [0, 0], [0, 1]]))
        (1 0
         0 3)
        """
        assert self.n == self.m
        assert space.n == self.n or space.m == 0
        assert space.rang() == space.m
        vectors = [space.find_inverse(self * space[:, j]) for j in range(space.m)]
        assert all(vectors[j] is not None for j in range(space.m))
        return Matrix([[vectors[j][i, 0] for j in range(space.m)] for i in range(space.m)])

    def find_base(self, n = None):
        """
        >>> Matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).find_base()
        (1 2 1
         4 5 0
         7 8 0)
        """
        if n is None:
            n = self.n
        reduced = self.append_right(Matrix.identity(n))
        j = 1
        while j <= reduced.m:
            if reduced[:, :j].rang() == j:
                j += 1
            else:
                reduced = reduced[:, :j - 1].append_right(reduced[:, j:])
        return reduced

    def eigenspace(self, eigenvalue):
        assert self.n == self.m
        return (self - Matrix.identity(self.n) * eigenvalue).null()

    def generalized_eigenspace(self, eigenvalue):
        assert self.n == self.m
        return ((self - Matrix.identity(self.n) * eigenvalue)**self.n).null()

    def nilpotent(self):
        return self.generalized_eigenspace(0).m == self.n

    def find_eigenvalue(self):
        """
        >>> Matrix([[1, 2, 0], [0, 1, 0], [0, 0, 3]]).find_eigenvalue()
        1
        """
        assert self.n == self.m
        for eigenvalue in range(-10, 11):
            if self.eigenspace(eigenvalue).m != 0:
                return eigenvalue
        raise RuntimeError("Unable to find an eigenvalue")

    def triangularize_1(self):
        """
        >>> (m := Matrix([[1, 2, 0], [0, 1, 0], [0, 0, 3]]).change_basis(Matrix([[3, 2, 5], [1, 4, 2], [5, 5, 1]])))
        (247/75 208/75 56/75
         -8/75  13/75  -34/75
         -14/15 4/15   23/15 )
        >>> (b := m.triangularize_1())
        (16  2  1
         1   0  0
         -10 -1 0)
        >>> m.change_basis(b)
        (3 6/25 -8/75
         0 1    2
         0 0    1    )
        """
        assert self.n == self.m
        if self.n == 0:
            return Matrix([])
        else:
            eigenvalue = self.find_eigenvalue()
            space = (self - Matrix.identity(self.n) * eigenvalue).range()
            triangularized = self.restrict(space).triangularize_1()
            return (space * triangularized).find_base(self.n)

    def triangularize_2(self):
        """
        >>> (m := Matrix([[1, 2, 0], [0, 1, 0], [0, 0, 3]]).change_basis(Matrix([[3, 2, 5], [1, 4, 2], [5, 5, 1]])))
        (247/75 208/75 56/75
         -8/75  13/75  -34/75
         -14/15 4/15   23/15 )
        >>> (b := m.triangularize_2())
        (2  -1 1
         -3 1  0
         5  0  0)
        >>> m.change_basis(b)
        (1 6/25 -14/75
         0 1    -2/3
         0 0    3     )
        """
        assert self.n == self.m
        if self.n == 0:
            return Matrix([])
        else:
            eigenvalue = self.find_eigenvalue()
            eigenvector = self.eigenspace(eigenvalue)[:, 0]
            space = eigenvector.find_base(self.n)
            triangularized = self.change_basis(space)[1:, 1:].triangularize_2()
            padded = Matrix([[1]]).append_bottom_right(triangularized)
            return space * padded

    def diagonal(self):
        assert self.n == self.m
        return [self[i, i] for i in range(self.n)]

    def eigenvalues(self):
        return self.change_basis(self.triangularize_2()).diagonal()

    def collected_eigenvalues(self):
        return list(Counter(self.eigenvalues()).items())

    def characteristic_polynomial(self):
        """
        >>> Matrix([[0, 0, -4], [1, 2, 2], [-1, 0, 0]]).characteristic_polynomial()
        8 * x^0 + -4 * x^1 + -2 * x^2 + 1 * x^3
        """
        assert self.n == self.m
        return prod((Polynomial([-x, 1]) for x in self.eigenvalues()), start=Polynomial([1]))

    def minimal_polynomial(self):
        """
        >>> Matrix([[0, 0, -4], [1, 2, 2], [-1, 0, 0]]).minimal_polynomial()
        -4 * x^0 + 0 * x^1 + 1 * x^2
        """
        assert self.n == self.m
        for n in range(self.n + 1):
            null = Matrix([list(self**i) for i in range(n)]).transpose().null()
            if null.m != 0:
                return Polynomial([null[i, 0] for i in range(n)])
        assert False

    def jordan_basis(self):
        """
        >>> (m := Matrix([[1, 2, 0], [0, 1, 0], [0, 0, 3]]).change_basis(Matrix([[3, 2, 5], [1, 4, 2], [5, 5, 1]])))
        (247/75 208/75 56/75
         -8/75  13/75  -34/75
         -14/15 4/15   23/15 )
        >>> (b := m.jordan_basis())
        (-2 25/6  -16
         3  -25/6 -1
         -5 0     10 )
        >>> m.change_basis(b)
        (1 1 0
         0 1 0
         0 0 3)
        >>> (m := Matrix([[5, 1, 0, 0, 0, 0], [0, 5, 1, 0, 0, 0], [0, 0, 5, 0, 0, 0], [0, 0, 0, 5, 1, 0], [0, 0, 0, 0, 5, 0], [0, 0, 0, 0, 0, 5]]).change_basis(Matrix([[3, 2, 5, 1, 6, 4], [1, 4, 2, 4, 1, 3], [5, 5, 1, 7, 9, 6], [4, 6, 8, 1, 2, 3], [6, 7, 9, 1, 2, 3], [8, 4, 2, 7, 4, 3]])))
        (-207/188  -1267/188 -1943/188  175/188   -87/188   -177/94
         -11/752   3071/752  671/752    -1843/752 -1119/752 -285/188
         4675/752  5561/752  11865/752  -309/752  311/752   429/188
         5741/752  6463/752  9631/752   3181/752  601/752   489/188
         3711/752  4109/752  8253/752   -2809/752 1715/752  -39/188
         -9211/752 -9833/752 -19049/752 5893/752  3537/752  889/188 )
        >>> (b := m.jordan_basis())
        (720   -2350/7  727795/567   -405/14 43705/756  1855975/18144
         -1380 0        -168025/378  60      -29875/504 -60775/432
         -60   0        -493735/1134 0       11075/1512 23375/1296
         -660  2820/7   0            180/7   0          4675/112
         -1020 -6110/7  0            585/14  0          -60775/672
         2860  23500/21 0            -830/7  0          116875/1008  )
        >>> m.change_basis(b)
        (5 1 0 0 0 0
         0 5 1 0 0 0
         0 0 5 0 0 0
         0 0 0 5 1 0
         0 0 0 0 5 0
         0 0 0 0 0 5)
        """
        assert self.n == self.m
        if self.nilpotent():
            def chains_to_matrix(chains):
                return reduce(Matrix.append_right, (vector for chain in chains for vector in chain), Matrix([]))
            def inner(self):
                if self.n == 0:
                    return []
                elif self.n == 1:
                    return [[Matrix([[1]])]]
                else:
                    space = self.range()
                    chains = inner(self.restrict(space))
                    for chain in chains:
                        for i in range(len(chain)):
                            chain[i] = space * chain[i]
                        chain.append(self.find_inverse(chain[-1]))
                    chains_matrix = chains_to_matrix(chains)
                    chains_matrix_base = chains_matrix.find_base()
                    chains_matrix_base_vectors = [chains_matrix_base[:, i] for i in range(chains_matrix.m, self.n)]
                    for chains_matrix_base_vector in chains_matrix_base_vectors:
                        chains.append([chains_matrix_base_vector - chains_matrix * Matrix([[0]]).append_bottom(chains_matrix.find_inverse(self * chains_matrix_base_vector)[:-1, :])])
                    return chains
            return chains_to_matrix(inner(self))
        else:
            result = Matrix([])
            for eigenvalue, _ in self.collected_eigenvalues():
                space = self.generalized_eigenspace(eigenvalue)
                result = result.append_right(space * ((self - Matrix.identity(self.n) * eigenvalue).restrict(space)).jordan_basis())
            return result

if __name__ == "__main__":
    import doctest
    doctest.testmod()
