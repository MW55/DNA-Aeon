import numpy
import typing

from norec4dna.distributions.Distribution import Distribution


class IdealSolitonDistribution(Distribution):
    def __init__(self, S: int = 5, seed: int = 0):  # S = Anzahl Bloecke
        super().__init__()
        self.rng = numpy.random
        self.rng.seed(seed)
        self.S: int = S
        self.pre_comp_dist: typing.List[float] = self.preCompute(S)

    def get_config_string(self) -> str:
        return "LT_IdealSoliton"

    def getNumber(self, seed: typing.Optional[int] = None) -> int:
        if seed is not None:
            self.set_seed(seed)
        return self.rng.choice(numpy.arange(1, self.S), p=self.pre_comp_dist)

    def preCompute(self, N: int) -> typing.List[float]:
        return self.normalize([1 / N] + [1 / (d * (d - 1)) for d in range(2, N)])

    def update_number_of_chunks(self, num_chunks):
        self.S = num_chunks
        self.pre_comp_dist = self.preCompute(self.S)


if __name__ == "__main__":
    x = IdealSolitonDistribution(20, 123)
    print(sum(x.pre_comp_dist))
    print([x.getNumber() for _ in range(1, 100)])
    print(x.get_distribution())
