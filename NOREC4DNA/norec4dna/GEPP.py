# Partially Based on jgcastro89 's Code https://gist.github.com/jgcastro89/49090cc69a499a129413597433b9baab
import numpy as np
import typing

from norec4dna.helper import xor_numpy

debug = False


def GEPP(a: np.array, b: np.array):
    return GEPP_intern(a, np.frombuffer(b, dtype="uint8"))


class GEPP_intern:
    """
    Gaussian elimination with partial pivoting.
    input: A is an n x n np matrix
           b is an n x 1 np array
    output: x is the solution of Ax=b
            with the entries permuted in
            accordance with the pivoting
            done by the algorithm
    post-condition: A and b have been modified.
    :return
    """

    def __init__(self, A: np.array, b: np.array):
        self.A: np.array = A  # input: A is an n x n np matrix
        self.b: np.array = b  # np.fromstring(b, dtype='uint8')  # b is an n x 1 np array
        self.packet_mapping: np.array = np.array([x + 1 for x in range(len(self.A))], dtype=np.uint)
        self.n: int = 0  # n is the length of A
        self.m: int = 0  # m is the width of A
        self.tmp_A: typing.List = []
        self.tmp_b: typing.List = []
        self._update_input()  # method that validates input
        self.result_mapping: np.array = np.zeros((np.int64(self.m), np.int64(1)), np.int32)

    def solve(self, partial=False) -> bool:
        if not (self.isPotentionallySolvable() or partial):
            return False
        try:
            self.insert_tmp()
            # since we need rows >= columns to perform partial pivoting, we fill with full-false rows..
            diff = len(self.A[0]) - len(self.A)
            if diff > 0 and partial:
                for _ in range(diff):
                    self.A = np.vstack((self.A, np.full((1, len(self.A[0])), False, bool)))
                    self.b = np.vstack((self.b, np.array([0 for x in range(len(self.b[0]))], dtype="uint8")))
                    self.packet_mapping = np.concatenate((self.packet_mapping, [len(self.A) + 1]))
                self._update_input()
            return self._elimination()
        except Exception as ex:
            print(f"GEPP: {ex}")
            return False

    def addRow(self, row: typing.Union[typing.Set[int], np.array], data: np.array, ):
        if self.n == 0:
            # first insert
            self.A = np.vstack((self.A, row))
            self.b = np.vstack((self.b, data))
            self.packet_mapping = np.vstack((self.packet_mapping, [len(self.A)]))
            self._update_input()
        self.tmp_A.append(row)
        self.tmp_b.append(data)
        self._update_input()

    def insert_tmp(self):
        if len(self.tmp_A) > 0:
            self.packet_mapping = np.concatenate(
                (self.packet_mapping, np.array([x for x in range(len(self.A) + 1, self.m + len(self.tmp_A) + 1)])))
            self.A = np.vstack((self.A, self.tmp_A))
            self.tmp_A.clear()
        if len(self.tmp_b) > 0:
            self.b = np.vstack((self.b, self.tmp_b))
            self.tmp_b.clear()

    def _update_input(self):
        self.n = len(self.A) + len(self.tmp_A)  # Number of Line
        self.m = len(self.A[0])  # Number of Columns

    def isPotentionallySolvable(self):
        return self.m <= self.n

    def _py_elimination(self) -> bool:
        """
        k represents the current pivot row. Since GE traverses the matrix in the
        upper right triangle, we also use k for indicating the k-th diagonal
        column index.
        :return
        """
        # Create Order Vector
        # self.order = np.array([[i] for i in range(0, self.n)])
        # Elimination
        for k in range(self.m):  # - 1):
            # Pivot
            maxindex = abs(self.A[k:, k]).argmax() + k
            if self.A[maxindex, k] == 0:
                continue
                # Matrix is singular - raise ValueError("Matrix is singular.")
            # Swap
            if maxindex != k:
                self.A[[k, maxindex]] = self.A[[maxindex, k]]
                self.b[[k, maxindex]] = self.b[[maxindex, k]]
                tmp = self.packet_mapping[maxindex]
                self.packet_mapping[maxindex] = self.packet_mapping[k]
                self.packet_mapping[k] = tmp
            # Eliminate top to bottom:
            for row in range(k + 1, self.n):
                if self.A[row, k:][0] == 1:
                    self.A[row, k:] = xor_numpy(self.A[row, k:], self.A[k, k:])
                    # Same for data:
                    self.b[row] = xor_numpy(self.b[row], self.b[k])

        # Eliminate bottom to top:
        for k in range(self.m - 1, -1, -1):
            for row in range(self.n - 1, -1, -1):
                if self.A[row, k] == 1 and row != k:
                    self.A[row, k:] = xor_numpy(self.A[row, k:], self.A[k, k:])
                    # Same for data:
                    self.b[row] = xor_numpy(self.b[row], self.b[k])

        self.result_mapping = self.generateResultMapping()
        return self.isSolved()

    try:
        from cdnarules import elimination  # just to trigger exception before running...

        def _elimination(self) -> bool:
            from cdnarules import elimination
            elimination(self.A, self.b, self.packet_mapping)
            self.result_mapping = self.generateResultMapping()
            return self.isSolved()
    except:
        print("Gaussian Elimination - C Module failed to load, falling back to slow mode")

        def _elimination(self) -> bool:
            return self._py_elimination()

    def generateResultMapping(self) -> np.array:
        """
        returns which row maps to which raw-data-Chunk. 
        """
        res_rows = np.full((self.m, 1), -1, int)
        for k in range(self.n):
            row = self.A[k]
            if np.sum(row) == 1:
                x = np.where(row == True)[0][0]
                res_rows[x] = k
        return res_rows

    def isSolved(self) -> bool:
        all_true = set(x for x in range(self.m))
        res = len(all_true) == len(set(x[0] for x in self.result_mapping))
        if debug and not res:
            print("Chunk(s) " + str(all_true - set(x[0] for x in self.result_mapping)) + " are missing. ")
        return res

    def getSolvedCount(self) -> int:
        return len(set(x[0] for x in self.result_mapping))


def main():
    A = np.array([[0, 1, 0, 1], [0, 0, 1, 0], [0, 1, 1, 0], [0, 1, 0, 0], [1, 1, 1, 0]], dtype=bool)

    b = np.array(
        [
            [np.fromstring("Hallo", dtype="uint8")],
            [np.fromstring("Welt!", dtype="uint8")],
            [np.fromstring("Test3", dtype="uint8")],
            [np.fromstring("Test1", dtype="uint8")],
            [np.fromstring("12345", dtype="uint8")],
        ]
    )

    gauss_elim_piv = GEPP(np.copy(A), np.copy(b))

    mapping = gauss_elim_piv.generateResultMapping()
    print(mapping)
    gauss_elim_piv.solve()
    gauss_elim_piv.solve()
    print("Is solved: " + str(gauss_elim_piv.isSolved()))
    print(gauss_elim_piv.A)
    print([chr(gauss_elim_piv.b[i]) for i in gauss_elim_piv.result_mapping])
    # print([chr(x) for x in gauss_elim_piv.b[2][0]])


if __name__ == "__main__":
    main()
