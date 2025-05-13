#include "dijkstra.h"
#include "fib.h"
#include "nlopt.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define FIB
#define MP

#define DATA_NCOLS 4320
#define DATA_NROWS 2160
#define DATA_CELLS (DATA_NCOLS * DATA_NROWS)

// #define ANL_NCOLS 456
// #define ANL_NROWS 240
#define ANL_NCOLS 2164
#define ANL_NROWS 864
#define ANL_CELLS (ANL_NCOLS * ANL_NROWS)

#define ANL_ROW_START 480
#define ANL_COL_START 2160

#define TRAIN_SIZE 530
#define TRAIN_SIZE_SQ (TRAIN_SIZE * TRAIN_SIZE)

#define NUM_COEFFS 5

#define ORIGIN_CELL 0
#define DEST_CELL 1000

typedef struct {
  int source;
  int target;
  double time;
} train_row;

typedef struct {
  grid_def *grid;
  train_row *train_data;
  double *graph;
  double *sq_res;
  int *iteration;
} nlopt_data;

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

void read_train_data(train_row *train_data, const char *path) {
  FILE *fp = fopen(path, "rb");
  int source;
  int target;
  double time;

  for (size_t i = 0; i < TRAIN_SIZE_SQ; i++) {
    fread(&source, sizeof(int), 1, fp);
    fread(&target, sizeof(int), 1, fp);
    fread(&time, sizeof(double), 1, fp);
    train_row row = {source, target, time};
    train_data[i] = row;
  }
  fclose(fp);
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

int id_data_to_anl(int orig_id) {
  int row = (orig_id / DATA_NCOLS) - ANL_ROW_START;
  int col = (orig_id % DATA_NCOLS) - ANL_COL_START;
  int new_id = row * ANL_NCOLS + col;
  assert(new_id >= 0);
  return new_id;
}

int id_anl_to_data(int orig_id) {
  int row = (orig_id / ANL_NCOLS) + ANL_ROW_START;
  int col = (orig_id % ANL_NCOLS) + ANL_COL_START;
  int new_id = row * DATA_NCOLS + col;
  assert(new_id >= 0);
  return new_id;
}

coord id_data_to_coord(int id) {
  int row = id / DATA_NCOLS;
  int col = id % DATA_NCOLS;

  double lat = 90.0 - row * (1.0 / 12.0) - (1.0 / 24.0);
  double lon = -180.0 + col * (1.0 / 12.0) + (1.0 / 24.0);
  coord pos = {lat, lon};
  return pos;
}

double obj_f(unsigned n, const double *x, double *grad, void *f_data) {
  nlopt_data *data = f_data;
  int iteration = *(data->iteration);
  printf("Starting iteration %d\n", iteration);
  grid_def *grid = data->grid;
  train_row *train_data = data->train_data;
  double *graph = data->graph;
  double *sq_res = data->sq_res;

  long double sum_of_squares = 0;

  coeffs_5 current_coeffs = {x[0], x[1], x[2], x[3], x[4]};

  initialize_graph_5_sq(graph, grid, &current_coeffs);

#pragma omp parallel shared(graph, train_data, sq_res)
  {
    node *A[24];
    fib_heap *Q = make_fib_heap(23, A);

    double *node_dist = malloc(ANL_CELLS * sizeof(double));
    double edge_distances[8];
    int potential_neighbors[8];

#pragma omp barrier
#pragma omp for
    for (size_t i = 0; i < TRAIN_SIZE; i++) {
      // printf("Solving iteration %d\n", (int)i);
      int source_id = id_data_to_anl(train_data[i * TRAIN_SIZE].source);
      dijkstra_solve_fib(graph, node_dist, Q, edge_distances,
                         potential_neighbors, ANL_CELLS, ANL_NCOLS, source_id);
      for (size_t j = 0; j < TRAIN_SIZE; j++) {
        train_row current_row = train_data[i * TRAIN_SIZE + j];
        double target_time = current_row.time;
        int target_id = id_data_to_anl(current_row.target);
        double model_time = node_dist[target_id];
        // printf("Target time = %f, Model time = %f\n", target_time,
        // model_time);
        sq_res[i * TRAIN_SIZE + j] = pow((target_time - model_time), 2);
      }
    }

    free(Q);
    free(node_dist);
  }

  for (size_t i = 0; i < TRAIN_SIZE_SQ; i++) {
    sum_of_squares += sq_res[i];
  }
  *(data->iteration) += 1;
  if ((iteration % 10) == 0) {
    for (size_t i = 0; i < NUM_COEFFS; i++) {
      printf("Coeffs[%d] = %f\n", (int)i, x[i]);
    }
    printf("Result = %f\n", (double)sum_of_squares / TRAIN_SIZE_SQ);
  }

  return (double)sum_of_squares / TRAIN_SIZE_SQ;
}

double obj_out(unsigned n, const double *x, double *grad, void *f_data,
               char *path) {
  nlopt_data *data = f_data;
  int iteration = *(data->iteration);
  printf("Starting iteration %d\n", iteration);
  grid_def *grid = data->grid;
  train_row *train_data = data->train_data;
  double *graph = data->graph;
  double *sq_res = data->sq_res;

  long double sum_of_squares = 0;

  coeffs_4 current_coeffs = {x[0], x[1], x[2], x[3]};

  initialize_graph_4(graph, grid, &current_coeffs);

  FILE *fp = fopen(path, "w");
  fprintf(fp, "SOURCE_ID,TARGET_ID,EXP_TIME,RES_TIME\n");

  node *A[24];
  fib_heap *Q = make_fib_heap(23, A);

  double *node_dist = malloc(ANL_CELLS * sizeof(double));
  double edge_distances[8];
  int potential_neighbors[8];

  for (size_t i = 0; i < TRAIN_SIZE; i++) {
    int source_id = train_data[i * TRAIN_SIZE].source;
    int source_id_local = id_data_to_anl(source_id);
    dijkstra_solve_fib(graph, node_dist, Q, edge_distances, potential_neighbors,
                       ANL_CELLS, ANL_NCOLS, source_id_local);
    for (size_t j = 0; j < TRAIN_SIZE; j++) {
      train_row current_row = train_data[i * TRAIN_SIZE + j];
      double target_time = current_row.time;
      int target_id = current_row.target;
      int target_id_local = id_data_to_anl(target_id);
      double model_time = node_dist[target_id_local];
      sq_res[i * TRAIN_SIZE + j] = pow((target_time - model_time), 2);
      fprintf(fp, "%d,%d,%f,%f\n", source_id, target_id, target_time,
              model_time);
    }
  }

  fclose(fp);
  free(Q);
  free(node_dist);

  for (size_t i = 0; i < TRAIN_SIZE_SQ; i++) {
    sum_of_squares += sq_res[i];
  }

  return (double)sum_of_squares / TRAIN_SIZE_SQ;
}

double obj_once(unsigned n, const double *x, double *grad, void *f_data,
                char *path) {
  nlopt_data *data = f_data;
  int iteration = *(data->iteration);
  printf("Starting iteration %d\n", iteration);
  grid_def *grid = data->grid;
  train_row *train_data = data->train_data;
  double *graph = data->graph;
  double *sq_res = data->sq_res;

  long double sum_of_squares = 0;

  coeffs_5 current_coeffs = {x[0], x[1], x[2], x[3], x[4]};

  initialize_graph_5_sq(graph, grid, &current_coeffs);

  node *A[24];
  fib_heap *Q = make_fib_heap(23, A);

  double *node_dist = malloc(ANL_CELLS * sizeof(double));
  int *node_prev = malloc(ANL_CELLS * sizeof(int));
  int *node_path = malloc(ANL_CELLS * sizeof(int));
  double edge_distances[8];
  int potential_neighbors[8];

  // printf("Solving iteration %d\n", (int)i);
  dijkstra_solve_fib_pred(graph, node_dist, node_prev, Q, edge_distances,
                          potential_neighbors, ANL_CELLS, ANL_NCOLS,
                          ORIGIN_CELL);

  FILE *fp = fopen(path, "w");
  fprintf(fp, "path\n");

  int prev_cell = DEST_CELL;
  while (prev_cell != ORIGIN_CELL) {
    int id_data = id_anl_to_data(prev_cell);
    coord prev_coord = id_data_to_coord(id_data);
    printf("%d -> (%f, %f)\n", id_data, prev_coord.lat, prev_coord.lon);
    fprintf(fp, "%d\n", id_data);
    prev_cell = node_prev[prev_cell];
  }

  fprintf(fp, "%d\n-1", id_anl_to_data(ORIGIN_CELL));
  fclose(fp);

  sum_of_squares += node_dist[DEST_CELL];
  free(Q);
  free(node_dist);
  free(node_prev);
  free(node_path);

  return sum_of_squares;
}

int main() {
  short *elev_full = calloc(DATA_CELLS, sizeof(short));
  short *elev_sub = calloc(ANL_CELLS, sizeof(short));

  short *terr_full = calloc(DATA_CELLS, sizeof(short));
  short *terr_sub = calloc(ANL_CELLS, sizeof(short));

  short *river_full = calloc(DATA_CELLS, sizeof(short));
  short *river_sub = calloc(ANL_CELLS, sizeof(short));

  read_data(elev_full, "./data/etopo1to5.bin");
  read_data(terr_full, "./data/landlake.bin");
  read_data(river_full, "./data/rivers.bin");

  grid_subset(elev_full, elev_sub, ANL_ROW_START, ANL_COL_START, ANL_NROWS,
              ANL_NCOLS, DATA_NCOLS);
  grid_subset(terr_full, terr_sub, ANL_ROW_START, ANL_COL_START, ANL_NROWS,
              ANL_NCOLS, DATA_NCOLS);
  grid_subset(river_full, river_sub, ANL_ROW_START, ANL_COL_START, ANL_NROWS,
              ANL_NCOLS, DATA_NCOLS);

  free(elev_full);
  free(terr_full);
  free(river_full);

  grid_def grid = {elev_sub,  terr_sub,      river_sub,    ANL_NCOLS,
                   ANL_NROWS, ANL_COL_START, ANL_ROW_START};

  train_row *train_data = calloc(TRAIN_SIZE_SQ, sizeof(train_row));
  read_train_data(train_data, "./data/train.bin");

  double *graph = calloc(ANL_CELLS * 8, sizeof(double));

  double *sq_res = calloc(TRAIN_SIZE_SQ, sizeof(double));

  int iteration = 0;

  nlopt_data data = {&grid, train_data, graph, sq_res, &iteration};

  double init_vals[NUM_COEFFS] = {0.288564, 0.0, 0.233367, 0.728704, 0.452882};
  double init_res = 11.674749;
#if 0
  nlopt_opt optimizer = nlopt_create(NLOPT_GN_MLSL_LDS, NUM_COEFFS);

  nlopt_opt local_optimizer = nlopt_create(NLOPT_LN_SBPLX, NUM_COEFFS);

  nlopt_result res_local =
      nlopt_set_local_optimizer(optimizer, local_optimizer);
  printf("Set local optimizer result code: %d\n", res_local);

  nlopt_result res_obj = nlopt_set_min_objective(optimizer, obj_f, &data);
  printf("Set objective result code: %d\n", res_obj);

  double lower_bounds[NUM_COEFFS] = {-10, -10, 0, 0, 0};
  nlopt_result res_lb = nlopt_set_lower_bounds(optimizer, lower_bounds);
  printf("Set lower bound result code: %d\n", res_lb);

  double upper_bounds[NUM_COEFFS] = {10, 10, 1, 1, 5};
  nlopt_result res_ub = nlopt_set_upper_bounds(optimizer, upper_bounds);
  printf("Set upper bound result code: %d\n", res_ub);

  nlopt_result res_maxeval = nlopt_set_maxeval(optimizer, 3000);
  printf("Set max eval param result code: %d\n", res_maxeval);

  nlopt_result res_opt = nlopt_optimize(optimizer, init_vals, &init_res);
  printf("Optimize result code: %d\n", res_opt);

  for (size_t i = 0; i < NUM_COEFFS; i++) {
    printf("Coeffs[%d] = %f\n", (int)i, init_vals[i]);
  }

  printf("Minimized result = %f\n", init_res);
#else
  // init_res = obj_f(NUM_COEFFS, init_vals, NULL, &data);
  // init_res = obj_out(NUM_COEFFS, init_vals, NULL, &data, "out.csv");
  init_res = obj_once(NUM_COEFFS, init_vals, NULL, &data, "path.csv");
  printf("Init result = %f\n", init_res);
#endif
  printf("Base Cost = %f\n", BASE_COST);
  return EXIT_SUCCESS;
}
