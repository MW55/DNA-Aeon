import numpy
import os
import struct
import glob
import typing
from math import ceil, floor
from reedsolo import RSCodec


def xor_mask(data, mask=0b11110000101011000011101001101010):
    return numpy.bitwise_xor(data, mask)


def get_file_size(file):
    return os.stat(file).st_size


class ReedSolomonEncoder:
    def __init__(self, file="logo.jpg", number_of_chunks=800, chunk_size=0, overhead=0.20):
        self.overhead: float = overhead
        self.chunk_size: int = chunk_size
        self.file: str = file
        self.insert_header: bool = True
        if self.chunk_size == 0:
            self.number_of_chunks: int = number_of_chunks
        else:
            self.set_no_chunks_from_chunk_size()
        file_size: int = get_file_size(self.file)
        self.chunk_size = ceil(1.0 * file_size / self.number_of_chunks)
        self.number_repair_symbols: int = int(floor((file_size // self.number_of_chunks) * self.overhead))
        self.rscodec: RSCodec = RSCodec(self.number_repair_symbols)

    def create_chunks(self) -> typing.List[bytes]:
        with open(self.file, "rb") as f:
            data = f.read()
        res = []
        for i in range(0, len(data), int(self.chunk_size)):
            encoded = self.rscodec.encode(struct.pack("<I" + str(len(data[i: i + self.chunk_size])) + "s", xor_mask(i),
                                                      data[i: i + self.chunk_size], ))
            res.append(struct.pack("<I" + str(len(encoded)) + "s", xor_mask(self.number_repair_symbols), encoded, ))
        self.number_of_chunks = len(res)
        return res

    def set_no_chunks_from_chunk_size(self):
        file_size: int = get_file_size(self.file)
        self.number_of_chunks = ceil(1.0 * file_size / self.chunk_size) + 1 if self.insert_header else 0
        print("number_of_chunks from chunk_size:" + str(self.number_of_chunks))

    def save_packets(self, split_to_multiple_files: bool, out_file: typing.Optional[str] = None):
        if not split_to_multiple_files:
            if out_file is None:
                out_file = self.file + ".RS"
            with open(out_file, "wb") as f:
                for packet in self.create_chunks():
                    f.write(packet)
        else:
            # Folder:
            if out_file is None:
                fulldir, filename = os.path.split(os.path.realpath(self.file))
                filename = "RS_" + filename
                out_file = os.path.join(fulldir, filename)

                if not out_file.endswith("/"):
                    files = glob.glob(out_file + "/*")
                else:
                    files = glob.glob(out_file + "*")
                for f in files:
                    os.remove(f)
            i = 0
            if not os.path.exists(out_file):
                os.makedirs(out_file)
            for packet in self.create_chunks():
                with open(out_file + "/" + str(i) + ".RS", "wb") as f:
                    f.write(packet)
                i += 1


class ReedSolomonDecoder:
    def __init__(self, file="logo.jpg"):
        self.decodeMap: typing.Dict[int, bytes] = {}
        self.file: str = file
        self.rscodec: typing.Optional[RSCodec] = None

    def decodeFolder(self):
        for file in os.listdir(self.file):
            if file.endswith(".RS"):
                with open(os.path.join(self.file, file), "rb") as f:
                    data = f.read()
                i, data = self.decode(data)
                self.decodeMap[i] = data
        with open(self.file + ".DECODED", "wb") as f:
            for key in sorted(self.decodeMap):
                f.write(self.decodeMap[key])

    def decode(self, text) -> typing.Tuple[int, bytes]:
        norepair_symbols: bytes = xor_mask(struct.unpack("<I", text[:4])[0])
        text = text[4:]
        if self.rscodec is None:
            self.rscodec = RSCodec(int(norepair_symbols))
        decoded = self.rscodec.decode(text)
        i: int = struct.unpack("<I", decoded[:4])[0]
        i: int = xor_mask(i)
        data: bytes = decoded[4:]
        return i, data

    @staticmethod
    def get_file_size(file) -> int:
        return os.stat(file).st_size


if __name__ == "__main__":
    rs = ReedSolomonEncoder("../.INFILES/logo.jpg", 500, 0, 0.20)
    rs.save_packets(True)