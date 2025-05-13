#include "fib.h"
#include <stdbool.h>
#include <stdio.h>

fib_heap *make_fib_heap(int max_degree, node **A) {
  fib_heap *empty_fib = malloc(sizeof(fib_heap));
  empty_fib->n = 0;
  empty_fib->min = NULL;
  empty_fib->max_degree = max_degree;
  empty_fib->A = A;
  return empty_fib;
}

void root_list_insert(fib_heap *H, node *new) {
  node *min = H->min;
  new->p = NULL;
  if (min->right == min) {
    min->left = new;
    new->right = min;
  } else {
    node *right_neighbor = min->right;
    new->right = right_neighbor;
    right_neighbor->left = new;
  }
  min->right = new;
  new->left = min;
}

void root_list_remove(fib_heap *H, node *old) {
  node *left_node = old->left;
  node *right_node = old->right;
  left_node->right = right_node;
  right_node->left = left_node;
}

void fib_insert(fib_heap *H, node *x) {
  x->degree = 0;
  x->p = NULL;
  x->child = NULL;
  x->mark = false;

  if (H->min == NULL) {
    H->min = x;
    x->right = x;
    x->left = x;
  } else {
    root_list_insert(H, x);
    if (x->key < H->min->key) {
      H->min = x;
    }
  }

  H->n += 1;
}

void child_list_insert(node *branch, node *root) {
  if (root->child == NULL) {
    root->child = branch;
    root->child->left = branch;
  } else {
    if (root->child == root->child->right) {
      root->child->left = branch;
      branch->right = root->child;
    } else {
      node *right_node = root->child->right;
      branch->right = right_node;
      right_node->left = branch;
    }
    branch->left = root->child;
  }
  root->child->right = branch;
  branch->p = root;
  root->degree += 1;
}

void child_list_remove(node *branch, node *parent) {
  assert(parent->degree > 0);

  if (parent->degree == 1) {
    parent->child = NULL;
  } else {
    node *left_node = branch->left;
    node *right_node = branch->right;
    left_node->right = right_node;
    right_node->left = left_node;
    if (parent->child == branch) {
      parent->child = right_node;
    }
  }

  parent->degree -= 1;
}

void merge_tree(fib_heap *H, node *root) {
  node *x = root;
  int d = x->degree;

  while (H->A[d] != NULL) {
    node *y = H->A[d];
    if (x->key > y->key) {
      node *temp = x;
      x = y;
      y = temp;
    }
    root_list_remove(H, y);
    child_list_insert(y, x);
    y->mark = false;
    H->A[d] = NULL;
    d += 1;
  }

  H->A[d] = x;
}

void fib_consolidate(fib_heap *H) {
  for (size_t i = 0; i <= H->max_degree; i++) {
    H->A[i] = NULL;
  }

  node *end_node = H->min->left;
  node *next_node = H->min;

  while (next_node != end_node) {
    node *current_node = next_node;
    next_node = current_node->right;
    merge_tree(H, current_node);
  }
  merge_tree(H, end_node);

  H->min = NULL;
  for (size_t i = 0; i <= H->max_degree; i++) {
    if (H->A[i] != NULL) {
      if (H->min == NULL) {
        H->min = H->A[i];
      } else {
        if (H->A[i]->key < H->min->key) {
          H->min = H->A[i];
        }
      }
    }
  }
}

node *fib_extract_min(fib_heap *H) {
  node *min = H->min;
  if (min != NULL) {
    if (min->child != NULL) {
      node *start = min->child;
      node *next = start->right;
      while (next != start) {
        node *current = next;
        next = next->right;
        root_list_insert(H, current);
      }
      root_list_insert(H, start);
    }
    if (min == min->right) {
      H->min = NULL;
    } else {
      root_list_remove(H, min);
      H->min = min->right;
      fib_consolidate(H);
    }
    H->n -= 1;
  }
  return min;
}

void cut_branch(fib_heap *H, node *branch, node *parent) {
  child_list_remove(branch, parent);
  root_list_insert(H, branch);
  branch->mark = false;
}

void cut_branch_cascade(fib_heap *H, node *branch) {
  node *parent = branch->p;

  if (parent != NULL) {
    if (branch->mark == false) {
      branch->mark = true;
    } else {
      cut_branch(H, branch, parent);
      cut_branch_cascade(H, parent);
    }
  }
}

void fib_decrease_key(fib_heap *H, node *x, double new_key) {
  assert(new_key < x->key);

  x->key = new_key;
  node *y = x->p;

  if ((y != NULL) && x->key < y->key) {
    cut_branch(H, x, y);
    cut_branch_cascade(H, y);
  }

  if (x->key < H->min->key) {
    H->min = x;
  }
}

void recursive_print(node *root) {
  if (root->child == NULL) {
    if (root->mark == true) {
      printf("(*%d)", root->id);
    } else {
      printf("(%d)", root->id);
    }
  } else {
    node *start = root->child;
    node *current = start;
    printf("([");
    while (current->right != start) {
      recursive_print(current);
      printf(",");
      current = current->right;
    }
    recursive_print(current);
    if (root->mark == true) {
      printf("]->*%d)", root->id);
    } else {
      printf("]->%d)", root->id);
    }
  }
}

void fib_print(fib_heap *H) {
  node *min = H->min;
  node *current = min;
  printf("----------\n");
  while (current->right != min) {
    recursive_print(current);
    printf("\n");
    current = current->right;
  }
  recursive_print(current);
  printf("\n");
  printf("----------\n");
}
