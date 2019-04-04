# Copyright 2019 Jeff Vierstra


GZIP_MAGIC_NUMBER="1f8b"

def is_gzip(filepath):
	with open(filepath) as f:
		if f.read(2).encode("hex")==GZIP_MAGIC_NUMBER:
			return True
		else:
			return False