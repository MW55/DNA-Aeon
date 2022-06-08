import numpy
import random
import typing
from math import ceil
from operator import add

from norec4dna.distributions.Distribution import Distribution


class RobustSolitonDistribution(Distribution):
    def __init__(self, S=10, K=8, delta=1.0, seed=0):
        # delta = failure probability # K = Grenze (position des hoehepunkts) , S = Anzahl der Bloecke
        super().__init__()
        self.rng: numpy.random = numpy.random
        self.rng.seed(seed)
        self.S: int = S
        self.K: int = K
        self.delta: float = delta
        self.pre_comp_dist: typing.List[float] = self.preCompute(S, K, delta)

    def get_config_string(self) -> str:
        return "LT_RobustSoliton_K=" + str(self.K) + "_delta=" + str(self.delta)

    @staticmethod
    def getGenerator(n: int) -> typing.Generator[int, typing.Any, typing.Any]:
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
    def robustSolitonDist(S: int, k: int, delta: float) -> typing.List[float]:
        R: float = S / k
        return ([1 / (d * k) for d in range(1, k)]
                + [numpy.log(R / delta) / k]
                + [0 for _ in range(k + 1, S)])

    def preCompute(self, S: int, k: int, delta: float) -> typing.List[float]:
        ideal = self.idealSolitonDist(S)
        robust = self.robustSolitonDist(S, k, delta)
        return self.normalize([i for i in map(add, ideal, robust)])

    def update_number_of_chunks(self, num_chunks: int):
        self.S = num_chunks
        self.pre_comp_dist = self.preCompute(self.S, self.K, self.delta)


if __name__ == "__main__":
    x = RobustSolitonDistribution(5, 2, 1.0, 123)
    print(sum(x.get_distribution()))
    print(x.get_distribution())
    print([x.getNumber() for _ in range(1, 100)])
