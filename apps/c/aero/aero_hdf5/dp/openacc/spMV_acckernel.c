//
// auto-generated by op2.py
//

//user function
int opDat0_spMV_stride_OP2CONSTANT;
int opDat0_spMV_stride_OP2HOST = -1;
int direct_spMV_stride_OP2CONSTANT;
int direct_spMV_stride_OP2HOST = -1;
//user function
//#pragma acc routine
inline void spMV_openacc( double **v, const double *K, const double **p) {

  v[0][0] += K[(0) * direct_spMV_stride_OP2CONSTANT] * p[0][0];
  v[0][0] += K[(1) * direct_spMV_stride_OP2CONSTANT] * p[1][0];
  v[1][0] += K[(1) * direct_spMV_stride_OP2CONSTANT] * p[0][0];
  v[0][0] += K[(2) * direct_spMV_stride_OP2CONSTANT] * p[2][0];
  v[2][0] += K[(2) * direct_spMV_stride_OP2CONSTANT] * p[0][0];
  v[0][0] += K[(3) * direct_spMV_stride_OP2CONSTANT] * p[3][0];
  v[3][0] += K[(3) * direct_spMV_stride_OP2CONSTANT] * p[0][0];
  v[1][0] += K[(4 + 1) * direct_spMV_stride_OP2CONSTANT] * p[1][0];
  v[1][0] += K[(4 + 2) * direct_spMV_stride_OP2CONSTANT] * p[2][0];
  v[2][0] += K[(4 + 2) * direct_spMV_stride_OP2CONSTANT] * p[1][0];
  v[1][0] += K[(4 + 3) * direct_spMV_stride_OP2CONSTANT] * p[3][0];
  v[3][0] += K[(4 + 3) * direct_spMV_stride_OP2CONSTANT] * p[1][0];
  v[2][0] += K[(8 + 2) * direct_spMV_stride_OP2CONSTANT] * p[2][0];
  v[2][0] += K[(8 + 3) * direct_spMV_stride_OP2CONSTANT] * p[3][0];
  v[3][0] += K[(8 + 3) * direct_spMV_stride_OP2CONSTANT] * p[2][0];
  v[3][0] += K[(15) * direct_spMV_stride_OP2CONSTANT] * p[3][0];
}

// host stub function
void op_par_loop_spMV(char const *name, op_set set,
  op_arg arg0,
  op_arg arg4,
  op_arg arg5){

  int nargs = 9;
  op_arg args[9];

  arg0.idx = 0;
  args[0] = arg0;
  for ( int v=1; v<4; v++ ){
    args[0 + v] = op_arg_dat(arg0.dat, v, arg0.map, 1, "double", OP_INC);
  }

  args[4] = arg4;
  arg5.idx = 0;
  args[5] = arg5;
  for ( int v=1; v<4; v++ ){
    args[5 + v] = op_arg_dat(arg5.dat, v, arg5.map, 1, "double", OP_READ);
  }


  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(3);
  op_timers_core(&cpu_t1, &wall_t1);
  OP_kernels[3].name      = name;
  OP_kernels[3].count    += 1;

  int  ninds   = 2;
  int  inds[9] = {0,0,0,0,-1,1,1,1,1};

  if (OP_diags>2) {
    printf(" kernel routine with indirection: spMV\n");
  }

  // get plan
  #ifdef OP_PART_SIZE_3
    int part_size = OP_PART_SIZE_3;
  #else
    int part_size = OP_part_size;
  #endif

  int set_size = op_mpi_halo_exchanges_cuda(set, nargs, args);


  int ncolors = 0;

  if (set_size >0) {

    if ((OP_kernels[3].count == 1) ||
        (opDat0_spMV_stride_OP2HOST != getSetSizeFromOpArg(&arg0))) {
      opDat0_spMV_stride_OP2HOST = getSetSizeFromOpArg(&arg0);
      opDat0_spMV_stride_OP2CONSTANT = opDat0_spMV_stride_OP2HOST;
    }
    if ((OP_kernels[3].count == 1) ||
        (direct_spMV_stride_OP2HOST != getSetSizeFromOpArg(&arg4))) {
      direct_spMV_stride_OP2HOST = getSetSizeFromOpArg(&arg4);
      direct_spMV_stride_OP2CONSTANT = direct_spMV_stride_OP2HOST;
    }

    //Set up typed device pointers for OpenACC
    int *map0 = arg0.map_data_d;

    double* data4 = (double*)arg4.data_d;
    double *data0 = (double *)arg0.data_d;
    double *data5 = (double *)arg5.data_d;

    op_plan *Plan = op_plan_get_stage(name,set,part_size,nargs,args,ninds,inds,OP_COLOR2);
    ncolors = Plan->ncolors;
    int *col_reord = Plan->col_reord;
    int set_size1 = set->size + set->exec_size;

    // execute plan
    for ( int col=0; col<Plan->ncolors; col++ ){
      if (col==1) {
        op_mpi_wait_all_cuda(nargs, args);
      }
      int start = Plan->col_offsets[0][col];
      int end = Plan->col_offsets[0][col+1];

      #pragma acc parallel loop independent deviceptr(col_reord,map0,data4,data0,data5)
      for ( int e=start; e<end; e++ ){
        int n = col_reord[e];
        int map0idx;
        int map1idx;
        int map2idx;
        int map3idx;
        map0idx = map0[n + set_size1 * 0];
        map1idx = map0[n + set_size1 * 1];
        map2idx = map0[n + set_size1 * 2];
        map3idx = map0[n + set_size1 * 3];

        double* arg0_vec[] = {
           &data0[1 * map0idx],
           &data0[1 * map1idx],
           &data0[1 * map2idx],
           &data0[1 * map3idx]};
        const double* arg5_vec[] = {
           &data5[1 * map0idx],
           &data5[1 * map1idx],
           &data5[1 * map2idx],
           &data5[1 * map3idx]};

        spMV_openacc(arg0_vec, &data4[n], arg5_vec);
      }

    }
    OP_kernels[3].transfer  += Plan->transfer;
    OP_kernels[3].transfer2 += Plan->transfer2;
  }

  if (set_size == 0 || set_size == set->core_size || ncolors == 1) {
    op_mpi_wait_all_cuda(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit_cuda(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[3].time     += wall_t2 - wall_t1;
}
