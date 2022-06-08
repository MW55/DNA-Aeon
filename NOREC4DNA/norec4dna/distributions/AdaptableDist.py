#!/usr/bin/python
import numpy
import typing

from norec4dna.distributions.Distribution import Distribution


class AdaptableDist(Distribution):
    """
    Allows Tools to directly set the desired Distribution
    E.g.: AdaptableDist
    """

    def __init__(self, dist: list, seed: int = 0):
        super().__init__()
        self.rng = numpy.random
        self.rng.seed(seed)
        self.S = len(dist) + 1
        self.dist = dist
        self.pre_comp_dist = self.normalize(dist)

    def getNumber(self, seed: typing.Optional[int] = None) -> int:
        if seed is not None:
            self.set_seed(seed)
        return self.rng.choice(numpy.arange(1, self.S), p=self.pre_comp_dist)

    def update_number_of_chunks(self, num_chunks: int):
        self.S = num_chunks
        self.pre_comp_dist = self.normalize(self.dist)
