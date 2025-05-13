#ifndef FIB_INCLUDED
#define FIB_INCLUDED
#include <assert.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>

typedef struct node {
  double key;
  int id;
  int degree;
  bool mark;
  struct node *p;
  struct node *child;
  struct node *left;
  struct node *right;
} node;

typedef struct {
  int n;
  node *min;
  int max_degree;
  node **A;
} fib_heap;

fib_heap *make_fib_heap(int max_degree, node **A);
void fib_insert(fib_heap *H, node *x);
node *fib_extract_min(fib_heap *H);
void fib_decrease_key(fib_heap *H, node *x, double new_key);
void fib_print(fib_heap *H);

#endif
