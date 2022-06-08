import numpy
import typing
from math import ceil, log

from norec4dna.distributions.Distribution import Distribution


class OnlineDistribution(Distribution):
    # Recommended parameters for large NumSourceBlocks: Epsilon=0.01, Quality=3
    # -> total_number_of_chunks = 1000 -> 3% overhead for recovery error rate of 1e-8
    def __init__(self, eps: float = 0.1, seed: int = 0):
        super().__init__()
        self.rng: numpy.random = numpy.random
        self.rng.seed(seed)
        self.eps: float = eps
        self.S: typing.Optional[int] = None
        self.pre_comp_dist: typing.List[float] = self.preCompute()

    def set_seed(self, seed: int):
        self.rng.seed(seed)

    def getNumber(self, seed: typing.Optional[int] = None) -> int:
        if seed is not None:
            self.set_seed(seed)
        return self.rng.choice(numpy.arange(1, self.S + 1), p=self.pre_comp_dist)

    def preCompute(self) -> typing.List[float]:
        s: int = ceil(log(self.eps * self.eps / 4) / log(1 - (self.eps / 2)))
        if self.S is None:
            self.S: int = s
        p1: float = 1 - ((1 + 1 / s) / (1 + self.eps))
        return self.normalize([p1] + [((1 - p1) * s) / ((s - 1) * i * (i - 1)) for i in range(2, s + 1)])

    def update_number_of_chunks(self, num_chunks: int):
        self.S = num_chunks
        self.pre_comp_dist = self.preCompute()


if __name__ == "__main__":
    x = OnlineDistribution(0.03, 123)
    print(sum(x.pre_comp_dist))
    print([x.getNumber() for _ in range(1, 100)])
    print(x.get_distribution())
