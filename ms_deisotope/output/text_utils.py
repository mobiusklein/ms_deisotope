import base64
import zlib

import numpy as np

from ms_deisotope.peak_set import Envelope, EnvelopePair

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
    encoded_string = base64.standard_b64encode(bytestring)
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


def envelopes_to_array(envelope_list, dtype=np.float32):
    collection = []
    for envelope in envelope_list:
        collection.append(0)
        collection.append(0)
        for pair in envelope:
            collection.extend(pair)
    return np.array(collection, dtype=np.float32)


def decode_envelopes(array):
    envelope_list = []
    current_envelope = []
    i = 0
    n = len(array)
    while i < n:
        a = array[i]
        b = array[i + 1]
        i += 2
        if a == 0 and b == 0:
            if current_envelope is not None:
                if current_envelope:
                    envelope_list.append(Envelope(current_envelope))
                current_envelope = []
        else:
            current_envelope.append(EnvelopePair(a, b))
    envelope_list.append(Envelope(current_envelope))
    return envelope_list


try:
    has_c = True
    _decode_envelopes = decode_envelopes
    _envelopes_to_array = envelopes_to_array
    from ms_deisotope._c.utils import decode_envelopes, envelopes_to_array
except ImportError:
    has_c = False
