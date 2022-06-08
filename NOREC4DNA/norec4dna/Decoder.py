import typing
import progressbar

from norec4dna.Packet import Packet


class Decoder:
    def __init__(self, file: typing.Optional[str] = None):
        self.file: typing.Optional[str] = file
        self.read_all_before_decode: bool = False
        self.isFolder: bool = False
        self.isZip: bool = False
        self.progress_bar = None

    @staticmethod
    def create_progress_bar(max_value):
        widgets = [progressbar.Percentage(), progressbar.Bar(), ' Correct: ', progressbar.Counter(), ', ',
                   progressbar.Variable('Corrupt'), ', ', progressbar.AdaptiveETA(), ' ', progressbar.Timer()]
        return progressbar.ProgressBar(max_value=max_value, widgets=widgets, max_error=False,
                                       redirect_stdout=True).start()

    @classmethod
    def pseudo_decoder(cls, number_of_chunks: typing.Optional[int] = None, read_all_before_decode: bool = False):
        def warn_pseudo():
            print("This method is not allowed while using pseudo Decoder")

        pseudo = cls(None)
        pseudo.decodeFile = warn_pseudo
        pseudo.decode = warn_pseudo
        pseudo.getNextValidPacket = warn_pseudo
        pseudo.saveDecodedFile = warn_pseudo
        pseudo.read_all_before_decode = read_all_before_decode
        if number_of_chunks is not None:
            pseudo.number_of_chunks = number_of_chunks
        pseudo.isPseudo = True
        return pseudo

    def input_new_packet(self, packet: Packet):
        pass  # implemented in subclasses

    def is_decoded(self):
        pass  # implemented in subclasses

    def solve(self):
        pass  # implemented in subclasses

    def decodeFolder(self, *args, **kwargs):
        pass  # implemented in subclasses

    def decodeFile(self, *args, **kwargs):
        pass  # implemented in subclasses

    def decodeZip(self, *args, **kwargs):
        pass  # implemented in subclasses

    def set_read_all_before_decode(self, do: bool):
        self.read_all_before_decode = do

    def decode(self, *args, **kwargs):
        if self.isFolder:
            return self.decodeFolder(*args, **kwargs)
        elif self.isZip:
            return self.decodeZip(*args, **kwargs)
        else:
            return self.decodeFile(*args, **kwargs)
