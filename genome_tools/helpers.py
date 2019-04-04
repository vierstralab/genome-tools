# Copyright 2019 Jeff Vierstra

import gzip

magic_dict = {
    b"\x1f\x8b\x08": "gz",
    b"\x42\x5a\x68": "bz2",
    b"\x50\x4b\x03\x04": "zip"
}

max_len = max(len(x) for x in magic_dict)

def get_file_type(filename):
    with open(filename, 'rb') as f:
        file_start = f.read(max_len)
        for magic, filetype in magic_dict.items():
            if file_start.startswith(magic):
                return filetype
        return None

def open_file(filename):
	file_type = get_file_type(filename)
	if file_type == "gz":
		return gzip.open(filename, mode='rt')
	else:
		return open(filename)