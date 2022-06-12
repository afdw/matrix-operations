import itertools
from fractions import Fraction
from math import gcd, lcm

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

    def __add__(self, other):
        assert(other.n == self.n and other.m == self.m)
        return Matrix([[self[i, j] + other[i, j] for j in range(self.m)] for i in range(self.n)])

    def __sub__(self, other):
        assert(other.n == self.n and other.m == self.m)
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
            assert(self.m == other.n)
            return Matrix([[sum(self[i, k] * other[k, j] for k in range(self.m)) for j in range(other.m)] for i in range(self.n)])

    def __truediv__(self, other):
        return Matrix([[self[i, j] / other for j in range(self.m)] for i in range(self.n)])

    def transpose(self):
        return Matrix([[self[i, j] for i in range(self.n)] for j in range(self.m)])

    def append_bottom(self, other):
        assert self.m == other.m
        return Matrix([[self[i, j] if i < self.n else other[i - self.n, j] for j in range(self.m)] for i in range(self.n + other.n)])

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

if __name__ == "__main__":
    import doctest
    doctest.testmod()
