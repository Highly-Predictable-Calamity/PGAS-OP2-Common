//
// auto-generated by op2.py
//

//user function
inline void res(const double *A, const float *u, float *du, const float *beta) {
  *du += (float)((*beta) * (*A) * (*u));
}
#ifdef VECTORIZE
//user function -- modified for vectorisation
inline void res_vec( const double *A, const float u[*][SIMD_VEC], float du[*][SIMD_VEC], const float *beta, int idx ) {
  du[0][idx]+= (float)((*beta) * (*A) * (u[0][idx]));

}
#endif

// host stub function
void op_par_loop_res(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3){

  int nargs = 4;
  op_arg args[4];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  //create aligned pointers for dats
  ALIGNED_double const double * __restrict__ ptr0 = (double *) arg0.data;
  __assume_aligned(ptr0,double_ALIGN);
  ALIGNED_float const float * __restrict__ ptr1 = (float *) arg1.data;
  __assume_aligned(ptr1,float_ALIGN);
  ALIGNED_float       float * __restrict__ ptr2 = (float *) arg2.data;
  __assume_aligned(ptr2,float_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(0);
  op_timers_core(&cpu_t1, &wall_t1);

  if (OP_diags>2) {
    printf(" kernel routine with indirection: res\n");
  }

  int exec_size = op_mpi_halo_exchanges(set, nargs, args);

  if (exec_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    for ( int n=0; n<(exec_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      float dat3[SIMD_VEC];
      for ( int i=0; i<SIMD_VEC; i++ ){
        dat3[i] = *((float*)arg3.data);
      }
      if (n+SIMD_VEC >= set->core_size) {
        op_mpi_wait_all(nargs, args);
      }
      ALIGNED_float float dat1[2][SIMD_VEC];
      ALIGNED_float float dat2[3][SIMD_VEC];
      #pragma omp simd simdlen(SIMD_VEC)
      for ( int i=0; i<SIMD_VEC; i++ ){
        int idx1_2 = 2 * arg1.map_data[(n+i) * arg1.map->dim + 1];

        dat1[0][i] = (ptr1)[idx1_2 + 0];
        dat1[1][i] = (ptr1)[idx1_2 + 1];

        dat2[0][i] = 0.0;
        dat2[1][i] = 0.0;
        dat2[2][i] = 0.0;

      }
      #pragma omp simd simdlen(SIMD_VEC)
      for ( int i=0; i<SIMD_VEC; i++ ){
        res_vec(
          &(ptr0)[3 * (n+i)],
          dat1,
          dat2,
          (float*)arg3.data,
          i);
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
        int idx2_3 = 3 * arg1.map_data[(n+i) * arg1.map->dim + 0];

        (ptr2)[idx2_3 + 0] += dat2[0][i];
        (ptr2)[idx2_3 + 1] += dat2[1][i];
        (ptr2)[idx2_3 + 2] += dat2[2][i];

      }
      for ( int i=0; i<SIMD_VEC; i++ ){
      }
    }

    //remainder
    for ( int n=(exec_size/SIMD_VEC)*SIMD_VEC; n<exec_size; n++ ){
    #else
    for ( int n=0; n<exec_size; n++ ){
    #endif
      if (n==set->core_size) {
        op_mpi_wait_all(nargs, args);
      }
      int map1idx = arg1.map_data[n * arg1.map->dim + 1];
      int map2idx = arg1.map_data[n * arg1.map->dim + 0];

      res(
        &(ptr0)[3 * n],
        &(ptr1)[2 * map1idx],
        &(ptr2)[3 * map2idx],
        (float*)arg3.data);
    }
  }

  if (exec_size == 0 || exec_size == set->core_size) {
    op_mpi_wait_all(nargs, args);
  }
  // combine reduction data
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[0].name      = name;
  OP_kernels[0].count    += 1;
  OP_kernels[0].time     += wall_t2 - wall_t1;
  OP_kernels[0].transfer += (float)set->size * arg1.size;
  OP_kernels[0].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[0].transfer += (float)set->size * arg0.size;
  OP_kernels[0].transfer += (float)set->size * arg3.size;
  OP_kernels[0].transfer += (float)set->size * arg1.map->dim * 4.0f;
}
