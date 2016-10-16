import base64
import zlib

import numpy as np


COMPRESSION_NONE = 'none'
COMPRESSION_ZLIB = 'zlib'


def encode_array(array, compression=COMPRESSION_NONE, dtype=np.float32):
    bytestring = np.asanyarray(array).astype(dtype).tobytes()
    if compression == COMPRESSION_NONE:
        bytestring = bytestring
    elif compression == COMPRESSION_ZLIB:
        bytestring = zlib.compress(bytestring)
    else:
        raise ValueError("Unknown compression: %s" % compression)
    encoded_string = base64.encodestring(bytestring)
    return encoded_string


def decode_array(bytestring, compression=COMPRESSION_NONE, dtype=np.float32):
    try:
        decoded_string = bytestring.encode("ascii")
    except AttributeError:
        decoded_string = bytestring
    decoded_string = base64.decodestring(decoded_string)
    if compression == COMPRESSION_NONE:
        decoded_string = decoded_string
    elif compression == COMPRESSION_ZLIB:
        decoded_string = zlib.decompress(decoded_string)
    else:
        raise ValueError("Unknown compression: %s" % compression)
    array = np.fromstring(decoded_string, dtype=dtype)
    return array
