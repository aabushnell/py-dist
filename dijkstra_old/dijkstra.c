#include "dijkstra.h"
#include <math.h>
#include <stdlib.h>

coord index_to_coord(size_t row, size_t col, bool centered,
                     const grid_def *grid) {
  double cell_size = 1.0 / 12.0;
  double offset = centered ? 0.5 : 0;
  double lat = 90.0 - cell_size * (grid->row_start + row + offset);
  double lon = -180.0 + cell_size * (grid->col_start + col + offset);
  coord c = {lat, lon};
  return c;
}

double make_rad(double degrees) {
  const double pi = acos(-1);
  return degrees * (pi / 180.0);
}

double hav(double angle) {
  double res = (1 - cos(angle)) / 2.0;
  return res;
}

double haversine_distance(size_t origin_row, size_t origin_col, size_t dest_row,
                          size_t dest_col, const grid_def *grid) {
  coord origin = index_to_coord(origin_row, origin_col, true, grid);
  coord dest = index_to_coord(dest_row, dest_col, true, grid);
  origin.lat = make_rad(origin.lat);
  dest.lat = make_rad(dest.lat);
  origin.lon = make_rad(origin.lon);
  dest.lon = make_rad(dest.lon);
  double res =
      12756.32 *
      asin(sqrt(hav(dest.lat - origin.lat) +
                (1 - hav(origin.lat - dest.lat) - hav(origin.lat + dest.lat)) *
                    hav(dest.lon - origin.lon)));
  return res;
}

double elevation_time_3(int elev_origin, int elev_dest,
                        const coeffs_3 *coeffs) {
  int avg = (int)((elev_dest + elev_origin) / 2);
  return avg * coeffs->elevation * 0.001; // convert m to km
}

double elevation_time_4(int elev_origin, int elev_dest,
                        const coeffs_4 *coeffs) {
  int avg = (int)((elev_dest + elev_origin) / 2);
  return avg * coeffs->elevation * 0.001; // convert m to km
}

double elevation_time_5_sq(int elev_origin, int elev_dest,
                           const coeffs_5 *coeffs) {
  int avg = (int)((elev_dest + elev_origin) / 2);
  return avg * coeffs->elevation * 0.001 +
         pow(avg, 2) * coeffs->elevation_sq * 0.001 * 0.001;
}

double elevation_time_7_up_down(int elev_origin, int elev_dest,
                                const coeffs_7 *coeffs) {
  int delta = elev_dest - elev_origin;
  if (delta >= 0) {
    return delta * coeffs->elev_up;
  } else if (delta > coeffs->elev_down_cutoff) {
    return delta * coeffs->elev_down;
  } else {
    return delta * coeffs->elev_down_steep;
  }
}

double elevation_time_full(int elev_origin, int elev_dest,
                           const time_coeffs *coeffs) {
  int delta = elev_dest - elev_origin;
  if (delta >= 0) {
    return delta * coeffs->up;
  } else if (delta < 0 && delta > coeffs->down_limit) {
    return delta * coeffs->down;
  } else {
    return delta * coeffs->down_steep;
  }
}

double travel_time_3(size_t origin_row, size_t origin_col, size_t dest_row,
                     size_t dest_col, const grid_def *grid,
                     const coeffs_3 *coeffs) {
  double distance =
      haversine_distance(origin_row, origin_col, dest_row, dest_col, grid);
  // printf("Distance = %f\n", distance);
  size_t origin_ind = origin_row * grid->ncols + origin_col;
  size_t dest_ind = dest_row * grid->ncols + dest_col;
  int origin_terrain = grid->terrain[origin_ind];
  int dest_terrain = grid->terrain[dest_ind];
  int origin_river = grid->river[origin_ind];
  int dest_river = grid->river[dest_ind];

  double time = 0.0;
  if (origin_terrain == 1 && dest_terrain == 1) {
    time = distance * BASE_COST + elevation_time_3(grid->elevation[origin_ind],
                                                   grid->elevation[dest_ind],
                                                   coeffs);
    if (origin_river == 1 && dest_river == 1) {
      time *= coeffs->river;
    }
  } else if (origin_terrain <= 0 && dest_terrain <= 0) {
    time = coeffs->sea * distance * BASE_COST;
  } else {
    if (origin_terrain <= 0) {
      time = distance * BASE_COST +
             elevation_time_3(0, grid->elevation[dest_ind], coeffs);
    } else {
      time = distance * BASE_COST +
             elevation_time_3(grid->elevation[origin_ind], 0, coeffs);
    }
  }
  return time;
}

double travel_time_4(size_t origin_row, size_t origin_col, size_t dest_row,
                     size_t dest_col, const grid_def *grid,
                     const coeffs_4 *coeffs) {
  double distance =
      haversine_distance(origin_row, origin_col, dest_row, dest_col, grid);
  // printf("Distance = %f\n", distance);
  size_t origin_ind = origin_row * grid->ncols + origin_col;
  size_t dest_ind = dest_row * grid->ncols + dest_col;
  int origin_terrain = grid->terrain[origin_ind];
  int dest_terrain = grid->terrain[dest_ind];
  int origin_river = grid->river[origin_ind];
  int dest_river = grid->river[dest_ind];

  double time = 0.0;
  if (origin_terrain == 1 && dest_terrain == 1) {
    time = distance * BASE_COST + elevation_time_4(grid->elevation[origin_ind],
                                                   grid->elevation[dest_ind],
                                                   coeffs);
    if (origin_river == 1 && dest_river == 1) {
      time *= coeffs->river;
    }
  } else if (origin_terrain <= 0 && dest_terrain <= 0) {
    time = coeffs->sea * distance * BASE_COST;
  } else {
    if (origin_terrain <= 0) {
      time = distance * BASE_COST +
             elevation_time_4(0, grid->elevation[dest_ind], coeffs) +
             coeffs->embark;
    } else {
      time = distance * BASE_COST +
             elevation_time_4(grid->elevation[origin_ind], 0, coeffs) +
             coeffs->embark;
    }
  }
  return time;
}

double travel_time_4_sq(size_t origin_row, size_t origin_col, size_t dest_row,
                        size_t dest_col, const grid_def *grid,
                        const coeffs_4 *coeffs) {
  double distance =
      haversine_distance(origin_row, origin_col, dest_row, dest_col, grid);
  // printf("Distance = %f\n", distance);
  size_t origin_ind = origin_row * grid->ncols + origin_col;
  size_t dest_ind = dest_row * grid->ncols + dest_col;
  int origin_terrain = grid->terrain[origin_ind];
  int dest_terrain = grid->terrain[dest_ind];
  int origin_river = grid->river[origin_ind];
  int dest_river = grid->river[dest_ind];

  double time = 0.0;
  if (origin_terrain == 1 && dest_terrain == 1) {
    time = distance * BASE_COST +
           pow(elevation_time_4(grid->elevation[origin_ind],
                                grid->elevation[dest_ind], coeffs),
               2);
    if (origin_river == 1 && dest_river == 1) {
      time *= coeffs->river;
    }
  } else if (origin_terrain <= 0 && dest_terrain <= 0) {
    time = coeffs->sea * distance * BASE_COST;
  } else {
    if (origin_terrain <= 0) {
      time = distance * BASE_COST +
             elevation_time_4(0, grid->elevation[dest_ind], coeffs) +
             coeffs->embark;
    } else {
      time = distance * BASE_COST +
             elevation_time_4(grid->elevation[origin_ind], 0, coeffs) +
             coeffs->embark;
    }
  }
  return time;
}

double travel_time_4_lim(size_t origin_row, size_t origin_col, size_t dest_row,
                         size_t dest_col, const grid_def *grid,
                         const coeffs_4 *coeffs) {
  double distance =
      haversine_distance(origin_row, origin_col, dest_row, dest_col, grid);
  // printf("Distance = %f\n", distance);
  size_t origin_ind = origin_row * grid->ncols + origin_col;
  size_t dest_ind = dest_row * grid->ncols + dest_col;
  int origin_terrain = grid->terrain[origin_ind];
  int dest_terrain = grid->terrain[dest_ind];
  int origin_river = grid->river[origin_ind];
  int dest_river = grid->river[dest_ind];
  int origin_elev = grid->elevation[origin_ind] / 1000;
  int dest_elev = grid->elevation[dest_ind] / 1000;

  double time = 0.0;
  if (origin_terrain == 1 && dest_terrain == 1) {
    if ((origin_elev < coeffs->elevation) && (dest_elev < coeffs->elevation)) {

      time = distance * BASE_COST;
      if (origin_river == 1 && dest_river == 1) {
        time *= coeffs->river;
      }
    } else {
      time = 1000 * 1000 * 1000;
    }
  } else if (origin_terrain <= 0 && dest_terrain <= 0) {
    time = coeffs->sea * distance * BASE_COST;
  } else {
    if (origin_terrain <= 0) {
      time = distance * BASE_COST + coeffs->embark;
    } else {
      time = distance * BASE_COST + coeffs->embark;
    }
  }
  return time;
}

double travel_time_5_sq(size_t origin_row, size_t origin_col, size_t dest_row,
                        size_t dest_col, const grid_def *grid,
                        const coeffs_5 *coeffs) {
  double distance =
      haversine_distance(origin_row, origin_col, dest_row, dest_col, grid);
  // printf("Distance = %f\n", distance);
  size_t origin_ind = origin_row * grid->ncols + origin_col;
  size_t dest_ind = dest_row * grid->ncols + dest_col;
  int origin_terrain = grid->terrain[origin_ind];
  int dest_terrain = grid->terrain[dest_ind];
  int origin_river = grid->river[origin_ind];
  int dest_river = grid->river[dest_ind];
  int origin_elev = grid->elevation[origin_ind];
  int dest_elev = grid->elevation[dest_ind];

  double time = 0.0;
  if (origin_terrain == 1 && dest_terrain == 1) {

    time = distance * BASE_COST +
           elevation_time_5_sq(origin_elev, dest_elev, coeffs);
    if (origin_river == 1 && dest_river == 1) {
      time *= coeffs->river;
    }
  } else if (origin_terrain <= 0 && dest_terrain <= 0) {
    time = coeffs->sea * distance * BASE_COST;
  } else {
    if (origin_terrain <= 0) {
      time = distance * BASE_COST + coeffs->embark +
             elevation_time_5_sq(0, dest_elev, coeffs);
    } else {
      time = distance * BASE_COST + coeffs->embark +
             elevation_time_5_sq(origin_elev, 0, coeffs);
    }
  }
  if (time >= 0) {
    return time;
  } else {
    return 0;
  }
}

double travel_time_7_up_down(size_t origin_row, size_t origin_col,
                             size_t dest_row, size_t dest_col,
                             const grid_def *grid, const coeffs_7 *coeffs) {
  double distance =
      haversine_distance(origin_row, origin_col, dest_row, dest_col, grid);
  // printf("Distance = %f\n", distance);
  size_t origin_ind = origin_row * grid->ncols + origin_col;
  size_t dest_ind = dest_row * grid->ncols + dest_col;
  int origin_terrain = grid->terrain[origin_ind];
  int dest_terrain = grid->terrain[dest_ind];
  int origin_river = grid->river[origin_ind];
  int dest_river = grid->river[dest_ind];
  int origin_elev = grid->elevation[origin_ind] / 1000;
  int dest_elev = grid->elevation[dest_ind] / 1000;

  double time = 0.0;
  if (origin_terrain == 1 && dest_terrain == 1) {

    time = distance * BASE_COST +
           elevation_time_7_up_down(origin_elev, dest_elev, coeffs);
    if (origin_river == 1 && dest_river == 1) {
      time *= coeffs->river;
    }
  } else if (origin_terrain <= 0 && dest_terrain <= 0) {
    time = coeffs->sea * distance * BASE_COST;
  } else {
    if (origin_terrain <= 0) {
      time = distance * BASE_COST + coeffs->embark +
             elevation_time_7_up_down(origin_elev, dest_elev, coeffs);
    } else {
      time = distance * BASE_COST + coeffs->embark +
             elevation_time_7_up_down(origin_elev, dest_elev, coeffs);
    }
  }
  if (time >= 0) {
    return time;
  } else {
    return 0;
  }
}

double travel_time_full(size_t origin_row, size_t origin_col, size_t dest_row,
                        size_t dest_col, const grid_def *grid,
                        const time_coeffs *coeffs) {
  double distance =
      haversine_distance(origin_row, origin_col, dest_row, dest_col, grid);
  // printf("Distance = %f\n", distance);
  size_t origin_ind = origin_row * grid->ncols + origin_col;
  size_t dest_ind = dest_row * grid->ncols + dest_col;
  int origin_terrain = grid->terrain[origin_ind];
  int dest_terrain = grid->terrain[dest_ind];
  int origin_river = grid->river[origin_ind];
  int dest_river = grid->river[dest_ind];

  double time = 0.0;
  if (origin_terrain == 1 && dest_terrain == 1) {
    time = distance + elevation_time_full(grid->elevation[origin_ind],
                                          grid->elevation[dest_ind], coeffs);
    if (origin_river == 1 && dest_river == 1) {
      time *= coeffs->river;
    }
  } else if (origin_terrain <= 0 && dest_terrain <= 0) {
    time = coeffs->sea * distance;
  } else {
    time = distance + coeffs->embark;
  }
  return time * BASE_COST;
}

void initialize_graph_3(double *graph, const grid_def *grid,
                        const coeffs_3 *coeffs) {
  size_t num_vertices = grid->ncols * grid->nrows;
  for (size_t v = 0; v < num_vertices; v++) {
    size_t row_ind = v / grid->ncols;
    size_t col_ind = v % grid->ncols;

    size_t offset = 0;
    for (int row_off = -1; row_off < 2; row_off++) {
      for (int col_off = -1; col_off < 2; col_off++) {
        if ((row_off == 0) && (col_off == 0)) {
          continue;
        }
        size_t new_row_ind = row_ind + row_off;
        size_t new_col_ind = col_ind + col_off;

        if (((new_row_ind >= 0) && (new_row_ind < grid->nrows)) &&
            ((new_col_ind >= 0) && (new_col_ind < grid->ncols))) {
          // printf("travel time between %d, %d and %d, %d\n", row_ind,
          // col_ind, new_row_ind, new_col_ind);
          double time = travel_time_3(row_ind, col_ind, new_row_ind,
                                      new_col_ind, grid, coeffs);
          graph[v * 8 + offset] = time;
        } else {
          graph[v * 8 + offset] = -1.0;
        }

        offset += 1;
      }
    }
  }
}

void initialize_graph_4(double *graph, const grid_def *grid,
                        const coeffs_4 *coeffs) {
  size_t num_vertices = grid->ncols * grid->nrows;
  for (size_t v = 0; v < num_vertices; v++) {
    size_t row_ind = v / grid->ncols;
    size_t col_ind = v % grid->ncols;

    size_t offset = 0;
    for (int row_off = -1; row_off < 2; row_off++) {
      for (int col_off = -1; col_off < 2; col_off++) {
        if ((row_off == 0) && (col_off == 0)) {
          continue;
        }
        size_t new_row_ind = row_ind + row_off;
        size_t new_col_ind = col_ind + col_off;

        if (((new_row_ind >= 0) && (new_row_ind < grid->nrows)) &&
            ((new_col_ind >= 0) && (new_col_ind < grid->ncols))) {
          // printf("travel time between %d, %d and %d, %d\n", row_ind,
          // col_ind, new_row_ind, new_col_ind);
          double time = travel_time_4(row_ind, col_ind, new_row_ind,
                                      new_col_ind, grid, coeffs);
          graph[v * 8 + offset] = time;
        } else {
          graph[v * 8 + offset] = -1.0;
        }

        offset += 1;
      }
    }
  }
}

void initialize_graph_4_sq(double *graph, const grid_def *grid,
                           const coeffs_4 *coeffs) {
  size_t num_vertices = grid->ncols * grid->nrows;
  for (size_t v = 0; v < num_vertices; v++) {
    size_t row_ind = v / grid->ncols;
    size_t col_ind = v % grid->ncols;

    size_t offset = 0;
    for (int row_off = -1; row_off < 2; row_off++) {
      for (int col_off = -1; col_off < 2; col_off++) {
        if ((row_off == 0) && (col_off == 0)) {
          continue;
        }
        size_t new_row_ind = row_ind + row_off;
        size_t new_col_ind = col_ind + col_off;

        if (((new_row_ind >= 0) && (new_row_ind < grid->nrows)) &&
            ((new_col_ind >= 0) && (new_col_ind < grid->ncols))) {
          // printf("travel time between %d, %d and %d, %d\n", row_ind,
          // col_ind, new_row_ind, new_col_ind);
          double time = travel_time_4_sq(row_ind, col_ind, new_row_ind,
                                         new_col_ind, grid, coeffs);
          graph[v * 8 + offset] = time;
        } else {
          graph[v * 8 + offset] = -1.0;
        }

        offset += 1;
      }
    }
  }
}

void initialize_graph_4_lim(double *graph, const grid_def *grid,
                            const coeffs_4 *coeffs) {
  size_t num_vertices = grid->ncols * grid->nrows;
  for (size_t v = 0; v < num_vertices; v++) {
    size_t row_ind = v / grid->ncols;
    size_t col_ind = v % grid->ncols;

    size_t offset = 0;
    for (int row_off = -1; row_off < 2; row_off++) {
      for (int col_off = -1; col_off < 2; col_off++) {
        if ((row_off == 0) && (col_off == 0)) {
          continue;
        }
        size_t new_row_ind = row_ind + row_off;
        size_t new_col_ind = col_ind + col_off;

        if (((new_row_ind >= 0) && (new_row_ind < grid->nrows)) &&
            ((new_col_ind >= 0) && (new_col_ind < grid->ncols))) {
          // printf("travel time between %d, %d and %d, %d\n", row_ind,
          // col_ind, new_row_ind, new_col_ind);
          double time = travel_time_4_lim(row_ind, col_ind, new_row_ind,
                                          new_col_ind, grid, coeffs);
          graph[v * 8 + offset] = time;
        } else {
          graph[v * 8 + offset] = -1.0;
        }

        offset += 1;
      }
    }
  }
}

void initialize_graph_5_sq(double *graph, const grid_def *grid,
                           const coeffs_5 *coeffs) {
  size_t num_vertices = grid->ncols * grid->nrows;
  for (size_t v = 0; v < num_vertices; v++) {
    size_t row_ind = v / grid->ncols;
    size_t col_ind = v % grid->ncols;

    size_t offset = 0;
    for (int row_off = -1; row_off < 2; row_off++) {
      for (int col_off = -1; col_off < 2; col_off++) {
        if ((row_off == 0) && (col_off == 0)) {
          continue;
        }
        size_t new_row_ind = row_ind + row_off;
        size_t new_col_ind = col_ind + col_off;

        if (((new_row_ind >= 0) && (new_row_ind < grid->nrows)) &&
            ((new_col_ind >= 0) && (new_col_ind < grid->ncols))) {
          // printf("travel time between %d, %d and %d, %d\n", row_ind,
          // col_ind, new_row_ind, new_col_ind);
          double time = travel_time_5_sq(row_ind, col_ind, new_row_ind,
                                         new_col_ind, grid, coeffs);
          graph[v * 8 + offset] = time;
        } else {
          graph[v * 8 + offset] = -1.0;
        }

        offset += 1;
      }
    }
  }
}

void initialize_graph_7_up_down(double *graph, const grid_def *grid,
                                const coeffs_7 *coeffs) {
  size_t num_vertices = grid->ncols * grid->nrows;
  for (size_t v = 0; v < num_vertices; v++) {
    size_t row_ind = v / grid->ncols;
    size_t col_ind = v % grid->ncols;

    size_t offset = 0;
    for (int row_off = -1; row_off < 2; row_off++) {
      for (int col_off = -1; col_off < 2; col_off++) {
        if ((row_off == 0) && (col_off == 0)) {
          continue;
        }
        size_t new_row_ind = row_ind + row_off;
        size_t new_col_ind = col_ind + col_off;

        if (((new_row_ind >= 0) && (new_row_ind < grid->nrows)) &&
            ((new_col_ind >= 0) && (new_col_ind < grid->ncols))) {
          // printf("travel time between %d, %d and %d, %d\n", row_ind,
          // col_ind, new_row_ind, new_col_ind);
          double time = travel_time_7_up_down(row_ind, col_ind, new_row_ind,
                                              new_col_ind, grid, coeffs);
          graph[v * 8 + offset] = time;
        } else {
          graph[v * 8 + offset] = -1.0;
        }

        offset += 1;
      }
    }
  }
}

void initialize_graph_full(double *graph, const grid_def *grid,
                           const time_coeffs *coeffs) {
  size_t num_vertices = grid->ncols * grid->nrows;
  for (size_t v = 0; v < num_vertices; v++) {
    size_t row_ind = v / grid->ncols;
    size_t col_ind = v % grid->ncols;

    size_t offset = 0;
    for (int row_off = -1; row_off < 2; row_off++) {
      for (int col_off = -1; col_off < 2; col_off++) {
        if ((row_off == 0) && (col_off == 0)) {
          continue;
        }
        size_t new_row_ind = row_ind + row_off;
        size_t new_col_ind = col_ind + col_off;

        if (((new_row_ind >= 0) && (new_row_ind < grid->nrows)) &&
            ((new_col_ind >= 0) && (new_col_ind < grid->ncols))) {
          // printf("travel time between %d, %d and %d, %d\n", row_ind,
          // col_ind, new_row_ind, new_col_ind);
          double time = travel_time_full(row_ind, col_ind, new_row_ind,
                                         new_col_ind, grid, coeffs);
          graph[v * 8 + offset] = time;
        } else {
          graph[v * 8 + offset] = -1.0;
        }

        offset += 1;
      }
    }
  }
}

size_t get_min_node(double *node_dist, bool *node_visited,
                    size_t num_vertices) {
  double min_dist = INFINITY;
  int min_node = -1;

  for (size_t v = 0; v < num_vertices; v++) {
    if ((node_dist[v] < min_dist) && (node_visited[v] == false)) {
      min_dist = node_dist[v];
      min_node = v;
    }
  }
  return min_node;
}

void validate_neighbors(size_t origin, const double *graph, size_t num_cols,
                        int *potential_neighbors, double *edge_distances) {
  size_t offset = 0;
  for (int row_off = -1; row_off < 2; row_off++) {
    for (int col_off = -1; col_off < 2; col_off++) {
      if (row_off == 0 && col_off == 0) {
        continue;
      }
      double edge_val = graph[origin * 8 + offset];
      if (edge_val < 0) {
        potential_neighbors[offset] = 0;
        edge_distances[offset] = -1.0;
      } else {
        potential_neighbors[offset] = row_off * num_cols + col_off;
        edge_distances[offset] = edge_val;
      }

      offset += 1;
    }
  }
}

void dijkstra_solve(const double *graph, double *node_dist, int *node_prev,
                    bool *node_visited, size_t num_vertices, size_t num_cols,
                    size_t origin) {
  // Initialize queue
  for (size_t v = 0; v < num_vertices; v++) {
    node_dist[v] = 1000 * 1000;
    node_prev[v] = -1;
    node_visited[v] = false;
  }

  node_dist[origin] = 0;

  // Traverse graph
  for (size_t i = 0; i < num_vertices; i++) {
    size_t u = get_min_node(node_dist, node_visited, num_vertices);
    node_visited[u] = true;
    int *potential_neighbors = calloc(8, sizeof(int));
    double *edge_distances = calloc(8, sizeof(double));
    validate_neighbors(u, graph, num_cols, potential_neighbors, edge_distances);
    for (size_t j = 0; j < 8; j++) {
      if (potential_neighbors[j] > 0) {
        size_t v = u + potential_neighbors[j];

        if (node_visited[v] == false) {
          double alt = node_dist[u] + edge_distances[j];
          if (alt < node_dist[v]) {
            node_dist[v] = alt;
            node_prev[v] = u;
          }
        }
      }
    }
    free(potential_neighbors);
    free(edge_distances);
  }
}

void dijkstra_solve_fib(const double *graph, double *node_dist, fib_heap *Q,
                        double *edge_distances, int *potential_neighbors,
                        size_t num_vertices, size_t num_cols, size_t origin) {
  node **node_addrs = malloc(num_vertices * sizeof(node *));
  node_dist[origin] = 0;

  for (size_t v = 0; v < num_vertices; v++) {
    if (v != origin) {
      node_dist[v] = 1000 * 1000 * 1000;
    }
    node *n_v = malloc(sizeof(node));
    // printf("Node allocated\n");
    n_v->key = node_dist[v];
    n_v->id = v;
    node_addrs[v] = n_v;
    fib_insert(Q, n_v);
    // printf("Node inserted\n");
  }
  while (Q->n != 0) {
    node *n_u = fib_extract_min(Q);
    // printf("Min is node %d, with dist %f\n", n_u->id, node_dist[n_u->id]);
    validate_neighbors(n_u->id, graph, num_cols, potential_neighbors,
                       edge_distances);

    for (size_t j = 0; j < 8; j++) {
      // printf("Edge Distance: %f\n", edge_distances[j]);
      if (edge_distances[j] > 0) {
        int v = n_u->id + potential_neighbors[j];
        // printf("Potential neighbor: %d\n", v);

        double alt = node_dist[n_u->id] + edge_distances[j];
        // printf("alt = %f\n", alt);
        if (alt < node_dist[v]) {
          // printf("Decrease Key!\n");
          node_dist[v] = alt;
          fib_decrease_key(Q, node_addrs[v], alt);
        }
      }
    }
    // break;
    free(n_u);
  }
  free(node_addrs);
}

void dijkstra_solve_fib_pred(const double *graph, double *node_dist,
                             int *node_prev, fib_heap *Q,
                             double *edge_distances, int *potential_neighbors,
                             size_t num_vertices, size_t num_cols,
                             size_t origin) {
  node **node_addrs = malloc(num_vertices * sizeof(node *));
  node_dist[origin] = 0;

  for (size_t v = 0; v < num_vertices; v++) {
    if (v != origin) {
      node_dist[v] = 1000 * 1000 * 1000;
      node_prev[v] = -1;
    }
    node *n_v = malloc(sizeof(node));
    n_v->key = node_dist[v];
    n_v->id = v;
    node_addrs[v] = n_v;
    fib_insert(Q, n_v);
  }

  while (Q->n != 0) {
    node *n_u = fib_extract_min(Q);
    validate_neighbors(n_u->id, graph, num_cols, potential_neighbors,
                       edge_distances);

    for (size_t j = 0; j < 8; j++) {
      if (edge_distances[j] > 0) {
        int v = n_u->id + potential_neighbors[j];

        double alt = node_dist[n_u->id] + edge_distances[j];
        if (alt < node_dist[v]) {
          node_dist[v] = alt;
          node_prev[v] = n_u->id;
          fib_decrease_key(Q, node_addrs[v], alt);
        }
      }
    }
    free(n_u);
  }
  free(node_addrs);
}
