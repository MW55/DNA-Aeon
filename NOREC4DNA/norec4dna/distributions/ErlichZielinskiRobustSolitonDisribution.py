import numpy
import random
import typing
from math import ceil
from operator import add

from norec4dna.distributions.Distribution import Distribution


class ErlichZielinskiRobustSolitonDistribution(Distribution):
    def __init__(self, k: int = 10, c: float = 0.025, delta: float = 0.001, seed: int = 0):
        # delta = failure probability # K = Grenze (position des hoehepunkts) , S = Anzahl der Bloecke
        super().__init__()
        self.rng = numpy.random
        self.set_seed(seed)
        self.c: float = c
        self.K: int = k
        self.S: int = k
        self.delta: float = delta
        self.pre_comp_dist: typing.List[float] = self.preCompute(k, c, delta)

    def update_number_of_chunks(self, num_chunks: int):
        self.K = num_chunks
        self.S = num_chunks
        self.pre_comp_dist = self.preCompute(self.K, self.c, self.delta)

    def get_config_string(self) -> str:
        return "LT_ErlichZielinskiRobustSoliton_K=" + str(self.K) + "_delta=" + str(self.delta) + "_c=" + str(self.c)

    @staticmethod
    def getGenerator(n: int):
        while True:
            num = random.random()
            res = int(ceil(1 / num))
            yield res if res <= n else 1

    def getNumber(self, seed: typing.Optional[int] = None) -> int:
        if seed is not None:
            self.set_seed(seed)
        return self.rng.choice(numpy.arange(1, self.S), p=self.pre_comp_dist)

    @staticmethod
    def idealSolitonDist(n: int) -> typing.List[float]:
        return [1 / n] + [1 / (d * (d - 1)) for d in range(2, n)]

    @staticmethod
    def robustSolitonDist(K: int, c: float, delta: float) -> typing.List[float]:
        s = c * numpy.sqrt(K) * numpy.power(numpy.log(1.0 * K / delta), 2)
        lim = int(round(1.0 * K / s))
        return ([s / (K * d) for d in range(1, lim)]
                + [(s * numpy.log(s / delta)) / K]
                + [0 for _ in range(lim + 1, K)]
                )

    def preCompute(self, k: int, c: float, delta: float) -> typing.List[float]:
        ideal = self.idealSolitonDist(k)
        robust = self.robustSolitonDist(k, c, delta)
        return self.normalize([i for i in map(add, ideal, robust)])


if __name__ == "__main__":
    x = ErlichZielinskiRobustSolitonDistribution(152, 0.025, 0.001, 123)
    print(sum(x.get_distribution()))
    print(x.get_distribution())
    print([x.getNumber() for _ in range(1, 100)])
