import json
import os
from pyteomics import xml, mzml, mzxml

try:
    basestring
except:
    basestring = (str, bytes)


def prepare_byte_index(byte_index):
    return byte_index.items()

def reconstruct_byte_index(mapping):
    index = xml.ByteEncodingOrderedDict()
    for key, value in mapping.items():
        index[key] = value
    return index

def _save_index(obj):
    file_name = obj.source_file
    if not isinstance(file_name, basestring):
        file_name = file_name.name
    index_name = file_name + '.idx'
    with open(index_name, 'wb') as handle:
        json.dump(prepare_byte_index(obj.index), handle)


def _load_index(obj):
    file_name = obj.source_file
    if not isinstance(file_name, basestring):
        file_name = file_name.name
    index_name = file_name + '.idx'
    with open(index_name, 'rb') as handle:
        mapping = json.load(handle)
    index = reconstruct_byte_index(mapping)
    return 

