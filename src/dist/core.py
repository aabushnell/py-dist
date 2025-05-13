import ctypes as C
import sys

import numpy as np

try:
    libcalc = C.CDLL('./lib/libcalc.so', mode=C.RTLD_GLOBAL)
except Exception as e:
    print('ERROR')
    print(e)
    sys.exit()

class MODEL_POINTER_STRUCT(C.Structure):
    _fields_ = [("edge_weight_neighbor_id", C.POINTER(C.c_int)),
                ("edge_weight_cost", C.POINTER(C.c_double)),
                ("node_data_node_id", C.POINTER(C.c_int)),
                ("node_data_count", C.POINTER(C.c_int)),
                ("node_data_expected", C.POINTER(C.c_int)),
                ("grid_data_elevation", C.POINTER(C.c_short)),
                ("grid_data_terrain", C.POINTER(C.c_short)),
                ("grid_data_rivers", C.POINTER(C.c_short)),
                ("grid_data_cell_weight", C.POINTER(C.c_double))]

alloc_model_raw = libcalc.alloc_model

alloc_model_raw.argtypes = [
    C.c_int,
    C.c_int,
    C.c_int
]

alloc_model_raw.restype = (
    MODEL_POINTER_STRUCT
)

class Model:

    def __init__(self, n_nodes: int, n_grid_cells: int, n_data_cells) -> None:
        self._n_nodes = n_nodes
        self._n_grid_cells = n_grid_cells
        self._n_data_cells = n_data_cells

        self._model_pointers = alloc_model_raw(n_nodes, n_grid_cells, n_data_cells)

        self.edge_weights = {
            'neighbor_ids' : np.ctypeslib.as_array(self._model_pointers.edge_weight_neighbor_id),
            'costs' : np.ctypeslib.as_array(self._model_pointers.edge_weight_cost)
        }

        self.node_data = {
            'node_ids' : np.ctypeslib.as_array(self._model_pointers.node_data_node_id),
            'counts' : np.ctypeslib.as_array(self._model_pointers.node_data_count),
            'expected' : np.ctypeslib.as_array(self._model_pointers.node_data_expected)
        }

        self.grid_data = {
            'elevations' : np.ctypeslib.as_array(self._model_pointers.grid_data_elevation),
            'terrain' : np.ctypeslib.as_array(self._model_pointers.grid_data_terrain),
            'rivers' : np.ctypeslib.as_array(self._model_pointers.grid_data_rivers),
            'cell_weights' : np.ctypeslib.as_array(self._model_pointers.grid_data_cell_weight)
        }






