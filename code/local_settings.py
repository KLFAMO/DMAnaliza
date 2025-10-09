project_path_str = ""
#project_path_str = "/home/pmorzynski/progs"

import pathlib as pa
if project_path_str:
    project_path = pa.Path(__file__).absolute().parents[2]
else:
    project_path = pa.Path(project_path_str)
