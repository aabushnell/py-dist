#ifndef DIJKSTRA_INCLUDED
#define DIJKSTRA_INCLUDED
#include "fib.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#define BASE_COST (1.0 / 30.0)

typedef struct {
  double up;
  double down;
  double down_steep;
  double down_limit;
  double sea;
  double river;
  double embark;
} time_coeffs;

typedef struct {
  double elevation;
  double sea;
  double river;
} coeffs_3;

typedef struct {
  double elevation;
  double sea;
  double river;
  double embark;
} coeffs_4;

typedef struct {
  double elevation;
  double elevation_sq;
  double sea;
  double river;
  double embark;
} coeffs_5;

typedef struct {
  double elev_up;
  double elev_down;
  double elev_down_steep;
  double elev_down_cutoff;
  double sea;
  double river;
  double embark;
} coeffs_7;

typedef struct {
  double lat;
  double lon;
} coord;

typedef struct {
  short *elevation;
  short *terrain;
  short *river;
  size_t ncols;
  size_t nrows;
  size_t col_start;
  size_t row_start;
} grid_def;

coord index_to_coord(size_t row, size_t col, bool centered,
                     const grid_def *grid);
double make_rad(double degrees);
double hav(double angle);

double haversine_distance(size_t origin_row, size_t origin_col, size_t dest_row,
                          size_t dest_col, const grid_def *grid);
double elevation_time_3(int elev_origin, int elev_dest, const coeffs_3 *coeffs);
double elevation_time_4(int elev_origin, int elev_dest, const coeffs_4 *coeffs);
double elevation_time_5_sq(int elev_origin, int elev_dest,
                           const coeffs_5 *coeffs);
double elevation_time_7_up_down(int elev_origin, int elev_dest,
                                const coeffs_7 *coeffs);
double elevation_time_full(int elev_origin, int elev_dest,
                           const time_coeffs *coeffs);

double travel_time_3(size_t origin_row, size_t origin_col, size_t dest_row,
                     size_t dest_col, const grid_def *grid,
                     const coeffs_3 *coeffs);
double travel_time_4(size_t origin_row, size_t origin_col, size_t dest_row,
                     size_t dest_col, const grid_def *grid,
                     const coeffs_4 *coeffs);
double travel_time_4_sq(size_t origin_row, size_t origin_col, size_t dest_row,
                        size_t dest_col, const grid_def *grid,
                        const coeffs_4 *coeffs);
double travel_time_4_lim(size_t origin_row, size_t origin_col, size_t dest_row,
                         size_t dest_col, const grid_def *grid,
                         const coeffs_4 *coeffs);
double travel_time_5_sq(size_t origin_row, size_t origin_col, size_t dest_row,
                        size_t dest_col, const grid_def *grid,
                        const coeffs_5 *coeffs);
double travel_time_7_up_down(size_t origin_row, size_t origin_col,
                             size_t dest_row, size_t dest_col,
                             const grid_def *grid, const coeffs_7 *coeffs);
double travel_time_full(size_t origin_row, size_t origin_col, size_t dest_row,
                        size_t dest_col, const grid_def *grid,
                        const time_coeffs *coeffs);

void initialize_graph_3(double *graph, const grid_def *grid,
                        const coeffs_3 *coeffs);
void initialize_graph_4(double *graph, const grid_def *grid,
                        const coeffs_4 *coeffs);
void initialize_graph_4_sq(double *graph, const grid_def *grid,
                           const coeffs_4 *coeffs);
void initialize_graph_4_lim(double *graph, const grid_def *grid,
                            const coeffs_4 *coeffs);
void initialize_graph_5_sq(double *graph, const grid_def *grid,
                           const coeffs_5 *coeffs);
void initialize_graph_7_up_down(double *graph, const grid_def *grid,
                                const coeffs_7 *coeffs);
void initialize_graph_full(double *graph, const grid_def *grid,
                           const time_coeffs *coeffs);

size_t get_min_node(double *node_dist, bool *node_visited, size_t num_vertices);
void validate_neighbors(size_t origin, const double *graph, size_t num_cols,
                        int *potential_neighbors, double *edge_distances);

void dijkstra_solve(const double *graph, double *node_dist, int *node_prev,
                    bool *node_visited, size_t num_vertices, size_t num_cols,
                    size_t origin);
void dijkstra_solve_fib(const double *graph, double *node_dist, fib_heap *Q,
                        double *edge_distances, int *potential_neighbors,
                        size_t num_vertices, size_t num_cols, size_t origin);
void dijkstra_solve_fib_pred(const double *graph, double *node_dist,
                             int *node_prev, fib_heap *Q,
                             double *edge_distances, int *potential_neighbors,
                             size_t num_vertices, size_t num_cols,
                             size_t origin);

#endif
