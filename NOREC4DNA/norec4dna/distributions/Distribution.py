import numpy as np
import typing


class Distribution(object):
    def __init__(self):
        self.rng = np.random
        self.pre_comp_dist: typing.List[float] = []
        self.S: typing.Optional[int] = None

    def get_config_string(self) -> str:
        return "Interface"

    @staticmethod
    def normalize(dist: typing.List[float]) -> typing.List[float]:
        return [float(i) / sum(dist) for i in dist]

    def get_distribution(self) -> typing.List[float]:
        return self.pre_comp_dist

    def get_size(self) -> int:
        return self.S

    def set_seed(self, seed: int):
        self.rng.seed(seed)

    def update_number_of_chunks(self, num_chunks: int):
        pass  # implemented in subclasses

    def getNumber(self, *args, **kwargs):
        pass
