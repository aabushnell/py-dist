#include "stdlib.h"

typedef struct {
  int *edge_weight_neighbor_id;
  double *edge_weight_cost;
  int *node_data_node_id;
  int *node_data_count;
  int *node_data_expected;
  short *grid_data_elevation;
  short *grid_data_terrain;
  short *grid_data_rivers;
  double *grid_data_cell_weight;
} model_pointer_struct;

model_pointer_struct alloc_model(int n_nodes, int n_grid_cells,
                                 int n_data_cells) {

  size_t data_size = n_nodes * 8 * 4         // edge_weight.neighbor_id (int)
                     + n_nodes * 8 * 8       // edge_weight.cost        (double)
                     + n_grid_cells * 1 * 4  // node_data.node_id       (int)
                     + n_grid_cells * 1 * 4  // node_data.count         (int)
                     + n_grid_cells * 1 * 4  // node_data.expected      (int)
                     + n_data_cells * 1 * 2  // grid_data.elevation     (short)
                     + n_data_cells * 1 * 2  // grid_data.terrain       (short)
                     + n_data_cells * 1 * 2  // grid_data.rivers        (short)
                     + n_data_cells * 1 * 8; // grid_data.cell_weight   (double)

  char *model_data = (char *)malloc(sizeof(char) * data_size);

  model_pointer_struct model_pointers;

  // clang-format off
  model_pointers.edge_weight_neighbor_id = (int *)(
    model_data + n_nodes * 0  + n_grid_cells * 0  + n_data_cells * 0);
  model_pointers.edge_weight_cost = (double *)(
    model_data + n_nodes * 32 + n_grid_cells * 0  + n_data_cells * 0);
  model_pointers.node_data_node_id = (int *)(
    model_data + n_nodes * 96 + n_grid_cells * 0  + n_data_cells * 0);
  model_pointers.node_data_count = (int *)(
    model_data + n_nodes * 96 + n_grid_cells * 4  + n_data_cells * 0);
  model_pointers.node_data_expected = (int *)(
    model_data + n_nodes * 96 + n_grid_cells * 8  + n_data_cells * 0);
  model_pointers.grid_data_elevation = (short *)(
    model_data + n_nodes * 96 + n_grid_cells * 12 + n_data_cells * 0);
  model_pointers.grid_data_terrain = (short *)(
    model_data + n_nodes * 96 + n_grid_cells * 12 + n_data_cells * 2);
  model_pointers.grid_data_rivers = (short *)(
    model_data + n_nodes * 96 + n_grid_cells * 12 + n_data_cells * 4);
  model_pointers.grid_data_cell_weight = (double *)(
    model_data + n_nodes * 96 + n_grid_cells * 12 + n_data_cells * 6);
  // clang-format on

  return model_pointers;
}
