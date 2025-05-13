#include "dijkstra.h"
#include "fib.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DATA_NCOLS 4320
#define DATA_NROWS 2160
#define DATA_CELLS (DATA_NCOLS * DATA_NROWS)

#define SEARCH_NCOLS 216
#define SEARCH_NROWS 216
#define SEARCH_CELLS (SEARCH_NCOLS * SEARCH_NROWS)

#define NODE_COUNT 1861
#define GRID_COUNT 3240

typedef struct {
  int node;
  double cost;
} edge;

typedef struct {
  int node_id;
  int count;
  int expected;
} node_data;

void read_data(short *grid, const char *path) {
  FILE *fp = fopen(path, "rb");
  const size_t ret_code = fread(grid, sizeof(grid[0]), DATA_CELLS, fp);
  if (ret_code != DATA_CELLS) {
    if (feof(fp)) {
      printf("ERROR: Unexpected end of file in input data\n");
    } else {
      perror("ERROR: Encountered reading input data");
    }
  }
  fclose(fp);
}

void read_weights(double *grid, const char *path) {
  FILE *fp = fopen(path, "rb");
  const size_t ret_code = fread(grid, sizeof(double), DATA_CELLS, fp);
  if (ret_code != DATA_CELLS) {
    if (feof(fp)) {
      printf("ERROR: Unexpected end of file in weights data\n");
    } else {
      perror("ERROR: Encountered reading weights data");
    }
  }
  fclose(fp);
}

int getfield(char *line, int num) {
  const char *tok;
  for (tok = strtok(line, ","); tok && *tok; tok = strtok(NULL, ",\n")) {
    if (!--num)
      return atoi(tok);
  }
  return -1;
}

int read_node_data(node_data *data, const char *path) {
  FILE *fp = fopen(path, "r");
  int count = 0;
  char line[100];
  fgets(line, 100, fp);
  printf("%s", line);
  while (fgets(line, 100, fp)) {
    char *tmp = strdup(line);
    char *tmp2 = strdup(line);
    data[count].node_id = getfield(tmp, 2);
    data[count].count = getfield(tmp2, 3);
    data[count].expected = getfield(line, 4);
    printf("%d, %d, %d\n", data[count].node_id, data[count].count,
           data[count].expected);
    free(tmp);
    free(tmp2);
    count++;
  }
  return count;
}

void write_edges(edge *edges, const char *path) {
  FILE *fp = fopen(path, "w");
  for (int i = 0; i < NODE_COUNT; i++) {
    for (int j = 0; j < 8; j++) {
      fprintf(fp, "%d,%d,%f\n", i, edges[i * 8 + j].node,
              edges[i * 8 + j].cost);
    }
  }
}

void grid_subset(short *orig, short *new, size_t row_start, size_t col_start,
                 size_t nrows, size_t ncols, size_t orig_cols) {
  for (size_t row = 0; row < nrows; row++) {
    for (size_t col = 0; col < ncols; col++) {
      new[row * ncols + col] =
          orig[(row + row_start) * orig_cols + col + col_start];
    }
  }
}

void double_subset(double *orig, double *new, size_t row_start,
                   size_t col_start, size_t nrows, size_t ncols,
                   size_t orig_cols) {
  for (size_t row = 0; row < nrows; row++) {
    for (size_t col = 0; col < ncols; col++) {
      new[row * ncols + col] =
          orig[(row + row_start) * orig_cols + col + col_start];
    }
  }
}

void apply_weights(double *results, double *distances, double *weights,
                   int starting_id, bool inland) {
  int size, center;
  if (inland) {
    size = 3;
    center = 1;
  } else {
    size = 9;
    center = 4;
  }
  int res_index = 0;

  for (int y = 0; y < size; y++) {
    for (int x = 0; x < size; x++) {
      if (y == center && x == center) {
        continue;
      }
      double sum = 0;
      for (int lil_y = y * 24; lil_y < (y + 1) * 24; lil_y++) {
        for (int lil_x = x * 24; lil_x < (x + 1) * 24; lil_x++) {
          int index = lil_y * size * 24 + lil_x;
          sum += distances[index] * weights[index];
        }
      }
      results[res_index++] += sum * weights[starting_id];
    }
  }
}

void place_sorted_edge(edge *edges, int node_id, int offset, edge candidate) {
  for (int i = offset; i < 8; i++) {
    double current_distance = edges[node_id * 8 + i].cost;
    if ((candidate.cost < current_distance && candidate.cost > 0) ||
        current_distance == 0) {
      if (current_distance > 0 && i < 7) {
        edge prev_edge = {edges[node_id * 8 + i].node, current_distance};
        place_sorted_edge(edges, node_id, i + 1, prev_edge);
      }
      edges[node_id * 8 + i] = candidate;
      break;
    }
  }
}

int choose_connections(int origin_id, double *weighted_distance, edge *edges,
                       node_data *nodes, bool inland) {

  int node_id = nodes[origin_id].node_id;
  int expected = nodes[origin_id].expected;
  int index = 0;
  int count = 0;
  int offset = inland ? -1 : 29;
  int offset_step_y = inland ? 0 : 6;

  // 00 01 02   03 04 05   06 07 08
  // 09 10 11   12 13 14   15 16 17
  // 18 19 20   21 22 23   24 25 26
  //          ------------
  // 27 28 29 | 30 31 32 | 33 34 35
  // 36 37 38 | 39 xx 40 | 41 42 43
  // 44 45 46 | 47 48 49 | 50 51 52
  //          ------------
  // 53 54 55   56 57 58   59 60 61
  // 62 63 64   65 66 67   68 69 70
  // 71 72 73   74 75 76   77 78 79

  for (int dy = -1; dy <= 1; dy++) {
    for (int dx = -1; dx <= 1; dx++) {
      // center (origin) tile
      if (dx == 0 && dy == 0) {
        continue;
      }
      offset++;
      int dest_id = origin_id + dy * 90 + dx;
      // skip out of bounds tiles
      if (dest_id < 0 || dest_id > 3239 ||
          (origin_id % 90 == 0 && dest_id % 90 == 89) ||
          (origin_id % 90 == 89 && dest_id % 90 == 0)) {
        continue;
      }
      int dest_node_id = nodes[dest_id].node_id;
      // skip sea tiles
      if (dest_node_id < 0) {
        continue;
      }
      edges[node_id * 8 + index] =
          (edge){dest_node_id, weighted_distance[offset]};
      count++;
      index++;
    }
    offset += offset_step_y;
  }

  while (index < 8 - expected + count) {
    edges[node_id * 8 + index] = (edge){-1, -1};
    index++;
  }
  if (index < 8) {
    offset = -1;
    for (int dy = -4; dy <= 4; dy++) {
      for (int dx = -4; dx <= 4; dx++) {
        // center (origin) tile
        if (dx == 0 && dy == 0) {
          continue;
        }
        offset++;
        // avoid double conting immediate neighbors
        if ((dy >= -1 && dy <= 1) && (dx >= -1 && dx <= 1)) {
          continue;
        }
        int dest_id = origin_id + dy * 90 + dx;
        // skip out of bounds tiles
        if (dest_id < 0 || dest_id > 3230 ||
            (origin_id % 90 < 4 && dest_id % 90 > 85) ||
            (origin_id % 90 > 85 && dest_id < 4)) {
          continue;
        }
        int dest_node_id = nodes[dest_id].node_id;
        // skip sea tiles
        if (dest_node_id < 0) {
          continue;
        }
        double distance = weighted_distance[offset];
        place_sorted_edge(edges, node_id, index,
                          (edge){dest_node_id, distance});
      }
    }
  }

  return count;
}

int main() {
  edge *edges = calloc(NODE_COUNT * 8, sizeof(edge));
  node_data *nodes = calloc(GRID_COUNT, sizeof(node_data));
  int nodes_counted = read_node_data(nodes, "./data/node_data.csv");
  if (nodes_counted != GRID_COUNT) {
    return EXIT_FAILURE;
  }

  size_t grid_data_size = DATA_CELLS * 3;
  short *grid_data = calloc(grid_data_size, sizeof(short));

  short *elevat_full = grid_data + 0 * DATA_CELLS;
  short *terran_full = grid_data + 1 * DATA_CELLS;
  short *rivers_full = grid_data + 2 * DATA_CELLS;

  read_data(elevat_full, "./data/etopo1to5.bin");
  read_data(terran_full, "./data/landlake.bin");
  read_data(rivers_full, "./data/rivers.bin");

  double *weight_full = calloc(DATA_CELLS, sizeof(double));
  read_weights(weight_full, "./data/weights.bin");

#pragma omp parallel shared(edges, nodes, elevat_full, terran_full,            \
                                rivers_full, weight_full)
  {
    size_t subs_data_size = SEARCH_CELLS * 3;
    short *subs_data = calloc(subs_data_size, sizeof(short));
    short *elevat_sub = subs_data + 0 * SEARCH_CELLS;
    short *terran_sub = subs_data + 1 * SEARCH_CELLS;
    short *rivers_sub = subs_data + 2 * SEARCH_CELLS;

    size_t graph_data_size = SEARCH_CELLS * 10;
    double *graph_data = calloc(graph_data_size, sizeof(double));
    double *graph = graph_data;
    double *node_dist = graph_data + 8 * SEARCH_CELLS;
    double *weight_sub = graph_data + 9 * SEARCH_CELLS;

    double edge_distances[8];
    int potential_neighbors[8];

    double weighted_distances[80];

    node *A[24];
    fib_heap *Q = make_fib_heap(23, A);

    coeffs_4 coeffs = {0.288564, 0.233367, 0.728704, 0.452882};
#pragma omp barrier
#pragma omp for
    for (int big_id = 0; big_id < GRID_COUNT; big_id++) {
      int tid = omp_get_thread_num();

      int big_y = big_id / 90;
      int big_x = big_id % 90;
      int search_start_y, search_start_x;
      int search_rows, search_cols, search_offset;
      bool inland = false;
      if (nodes[big_id].node_id < 0) {
        printf("THREAD %d: Sea Node at id %d\n", tid, big_id);
        continue;
      } else if (nodes[big_id].count == nodes[big_id].expected) {
        printf("THREAD %d: Inland Node at id %d\n", tid, big_id);
        inland = true;
        search_start_y = 216 + ((big_y - 1) * 24);
        search_start_x = 1848 + ((big_x - 1) * 24);
        search_rows = 3 * 24;
        search_cols = 3 * 24;
        search_offset = 24 * search_cols + 24;
      } else {
        printf("THREAD %d: Coastal Node at id %d\n", tid, big_id);
        search_start_y = 216 + ((big_y - 4) * 24);
        search_start_x = 1848 + ((big_x - 4) * 24);
        search_rows = 9 * 24;
        search_cols = 9 * 24;
        search_offset = 4 * 24 * search_cols + 4 * 24;
      }

      grid_subset(elevat_full, elevat_sub, search_start_y, search_start_x,
                  search_rows, search_cols, DATA_NCOLS);
      grid_subset(terran_full, terran_sub, search_start_y, search_start_x,
                  search_rows, search_cols, DATA_NCOLS);
      grid_subset(rivers_full, rivers_sub, search_start_y, search_start_x,
                  search_rows, search_cols, DATA_NCOLS);

      double_subset(weight_full, weight_sub, search_start_y, search_start_x,
                    search_rows, search_cols, DATA_NCOLS);

      grid_def grid = {elevat_sub,  terran_sub,     rivers_sub,    search_cols,
                       search_rows, search_start_x, search_start_y};

      initialize_graph_4(graph, &grid, &coeffs);

      for (int i = 0; i < 80; i++) {
        weighted_distances[i] = 0;
      }
      for (int lil_y = 0; lil_y < 24; lil_y++) {
        for (int lil_x = 0; lil_x < 24; lil_x++) {
          int lil_id = search_offset + lil_y * search_cols + lil_x;
          dijkstra_solve_fib(graph, node_dist, Q, edge_distances,
                             potential_neighbors, search_rows * search_cols,
                             search_cols, lil_id);
          apply_weights(weighted_distances, node_dist, weight_sub, lil_id,
                        inland);
        }
      }
      choose_connections(big_id, weighted_distances, edges, nodes, inland);
      printf("THREAD %d: Found neighbors at grid %d\n", tid, big_id);
    }
  }

  write_edges(edges, "./graph_edges.csv");

  return EXIT_SUCCESS;
}
