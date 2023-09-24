#progspath_str = "/home/pmorzynski/progs"
progspath_str = ""

import pathlib as pa
if progspath_str == "":
    progspath = pa.Path(__file__).absolute().parents[2]
else:
    progspath = pa.Path(progspath_str)