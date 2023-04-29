#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// global constants

double gam, gm1, cfl, eps, mach, alpha, qinf[4];

//
// OP header file
//

#include "op_seq.h"

//
// kernel routines for parallel loops
//

#include "do_nothing.h"
// main program

int main(int argc, char **argv) {
  // OP initialisation
  op_init(argc, argv, 2);

  int *becell, *ecell, *bound, *bedge, *edge, *cell;
  double *x, *q, *qold, *adt, *res;

  int nnode, ncell, nedge, nbedge, niter;
  double rms;

  // timer
  double cpu_t1, cpu_t2, wall_t1, wall_t2;

  // read in grid

  op_printf("reading in grid \n");

  FILE *fp;
  if ((fp = fopen("./new_grid_2880000.dat", "r")) == NULL) {
    op_printf("can't open file new_grid.dat\n");
    exit(-1);
  }

  if (fscanf(fp, "%d %d %d %d \n", &nnode, &ncell, &nedge, &nbedge) != 4) {
    op_printf("error reading from new_grid.dat\n");
    exit(-1);
  }

  cell = (int *)malloc(4 * ncell * sizeof(int));
  edge = (int *)malloc(2 * nedge * sizeof(int));
  ecell = (int *)malloc(2 * nedge * sizeof(int));
  bedge = (int *)malloc(2 * nbedge * sizeof(int));
  becell = (int *)malloc(nbedge * sizeof(int));
  bound = (int *)malloc(nbedge * sizeof(int));

  x = (double *)malloc(2 * nnode * sizeof(double));
  q = (double *)malloc(4 * ncell * sizeof(double));
  qold = (double *)malloc(4 * ncell * sizeof(double));
  res = (double *)malloc(4 * ncell * sizeof(double));
  adt = (double *)malloc(ncell * sizeof(double));

  for (int n = 0; n < nnode; n++) {
    if (fscanf(fp, "%lf %lf \n", &x[2 * n], &x[2 * n + 1]) != 2) {
      op_printf("error reading from new_grid.dat\n");
      exit(-1);
    }
  }

  for (int n = 0; n < ncell; n++) {
    if (fscanf(fp, "%d %d %d %d \n", &cell[4 * n], &cell[4 * n + 1],
               &cell[4 * n + 2], &cell[4 * n + 3]) != 4) {
      op_printf("error reading from new_grid.dat\n");
      exit(-1);
    }
  }

  for (int n = 0; n < nedge; n++) {
    if (fscanf(fp, "%d %d %d %d \n", &edge[2 * n], &edge[2 * n + 1],
               &ecell[2 * n], &ecell[2 * n + 1]) != 4) {
      op_printf("error reading from new_grid.dat\n");
      exit(-1);
    }
  }

  for (int n = 0; n < nbedge; n++) {
    if (fscanf(fp, "%d %d %d %d \n", &bedge[2 * n], &bedge[2 * n + 1],
               &becell[n], &bound[n]) != 4) {
      op_printf("error reading from new_grid.dat\n");
      exit(-1);
    }
  }

  fclose(fp);

  // set constants and initialise flow field and residual

  op_printf("initialising flow field \n");

  gam = 1.4f;
  gm1 = gam - 1.0f;
  cfl = 0.9f;
  eps = 0.05f;

  double mach = 0.4f;
  double alpha = 3.0f * atan(1.0f) / 45.0f;
  double p = 1.0f;
  double r = 1.0f;
  double u = sqrt(gam * p / r) * mach;
  double e = p / (r * gm1) + 0.5f * u * u;

  qinf[0] = r;
  qinf[1] = r * u;
  qinf[2] = 0.0f;
  qinf[3] = r * e;

  for (int n = 0; n < ncell; n++) {
    for (int m = 0; m < 4; m++) {
      q[4 * n + m] = qinf[m];
      res[4 * n + m] = 0.0f;
    }
  }

  // declare sets, pointers, datasets and global constants

  op_set nodes = op_decl_set(nnode, "nodes");
  op_set edges = op_decl_set(nedge, "edges");
  op_set bedges = op_decl_set(nbedge, "bedges");
  op_set cells = op_decl_set(ncell, "cells");

  op_map pedge = op_decl_map(edges, nodes, 2, edge, "pedge");
  free(edge);
  op_map pecell = op_decl_map(edges, cells, 2, ecell, "pecell");
  free(ecell);
  op_map pbedge = op_decl_map(bedges, nodes, 2, bedge, "pbedge");
  free(bedge);
  op_map pbecell = op_decl_map(bedges, cells, 1, becell, "pbecell");
  free(becell);
  op_map pcell = op_decl_map(cells, nodes, 4, cell, "pcell");
  free(cell);

  op_dat p_bound = op_decl_dat(bedges, 1, "int", bound, "p_bound");
  free(bound);
  op_dat p_x = op_decl_dat(nodes, 2, "double", x, "p_x");
  free(x);
  op_dat p_q = op_decl_dat(cells, 4, "double", q, "p_q");
  free(q);
  op_dat p_qold = op_decl_dat(cells, 4, "double", qold, "p_qold");
  free(qold);
  op_dat p_adt = op_decl_dat(cells, 1, "double", adt, "p_adt");
  free(adt);
  op_dat p_res = op_decl_dat(cells, 4, "double", res, "p_res");
  free(res);

  op_decl_const(1, "double", &gam);
  op_decl_const(1, "double", &gm1);
  op_decl_const(1, "double", &cfl);
  op_decl_const(1, "double", &eps);
  op_decl_const(1, "double", &mach);
  op_decl_const(1, "double", &alpha);
  op_decl_const(4, "double", qinf);

  printf("starting time-marching loop \n");
  fflush(stdout);

  op_diagnostic_output();

  printf("starting time-marching loop \n");
  fflush(stdout);

  // initialise timers for total execution wall time
  op_timers(&cpu_t1, &wall_t1);

  // main time-marching loop

  niter = 100000;

  printf("starting time-marching loop \n");
  fflush(stdout);

  for (int iter = 1; iter <= niter; iter++) {

    // Do nothing kernels to test performance of exchange on GPI
    op_par_loop(do_nothing, "do_nothing", bedges,
                  op_arg_dat(p_x, 0, pbedge, 2, "double", OP_READ),
                  op_arg_dat(p_x, 1, pbedge, 2, "double", OP_READ),
                  op_arg_dat(p_q, 0, pbecell, 4, "double", OP_READ),
                  op_arg_dat(p_adt, 0, pbecell, 1, "double", OP_READ),
                  op_arg_dat(p_res, 0, pbecell, 4, "double", OP_INC),
                  op_arg_dat(p_bound, -1, OP_ID, 1, "int", OP_READ));



    
  }

  op_timers(&cpu_t2, &wall_t2);

  // output the result dat array to files - no, I don't think so
//   op_print_dat_to_txtfile(p_q, "out_grid_seq.dat"); // ASCI
//   op_print_dat_to_binfile(p_q, "out_grid_seq.bin"); // Binary

  // write given op_dat's indicated segment of data to a memory block in the
  // order it was originally
  // arranged (i.e. before partitioning and reordering)
//   double *q_part = (double *)op_malloc(sizeof(double) * op_get_size(cells) * 4);
//   op_fetch_data_idx(p_q, q_part, 0, op_get_size(cells) - 1);
//   free(q_part);

  op_timing_output();
  op_printf("Max total runtime = %f\n", wall_t2 - wall_t1);

  op_exit();
}
