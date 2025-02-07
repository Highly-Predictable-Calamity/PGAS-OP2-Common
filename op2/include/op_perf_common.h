#include <op_lib_c.h>
#include <op_lib_mpi.h>
/*******************************************************************************
* Data structure to hold communications of an op_dat
*******************************************************************************/
#define NAMESIZE 20
typedef struct {
  // name of this op_dat
  char name[NAMESIZE];
  // size of this op_dat
  int size;
  // index of this op_dat
  int index;
  // total number of times this op_dat was halo exported
  int count;
  // total number of bytes halo exported for this op_dat in this kernel
  int bytes;
} op_dat_comm_info_core;

typedef op_dat_comm_info_core *op_dat_comm_info;

typedef struct{
    char name[NAMESIZE];
    double time;
    int count;
} kernel_breakdown;

/*******************************************************************************
* Data Type to hold performance measures for communication
*******************************************************************************/

typedef struct {

  UT_hash_handle hh; // with this variable uthash makes this structure hashable

  // name of kernel
  char name[NAMESIZE];
  // total time spent in this kernel (compute + comm - overlap)
  double time;
  // number of times this kernel is called
  int count;
  // number of op_dat indices in this kernel
  int num_indices;
  // array to hold all the op_dat_mpi_comm_info structs for this kernel
  op_dat_comm_info *comm_info;
  // capacity of comm_info array
  int cap;
  // array to hold extra timing breakdowns
  kernel_breakdown* breakdown;
  // capacity of the breakdown array
  int breakdown_cap;
} op_comm_kernel;

extern op_comm_kernel *op_comm_kernel_tab;

op_comm_kernel *get_kernel_entry(const char *name);
void *op_comm_perf_time(const char *name, double time);
kernel_breakdown *get_breakdown(op_comm_kernel *kernel_entry, const char *breakdown_name);
void *op_comm_perf_time_breakdown(const char *kernel_name, const char *breakdown_name, double time);
#ifdef COMM_PERF
int search_op_comm_kernel(op_dat dat, op_comm_kernel *kernal, int num_indices);
void op__perf_comm(void *k_i, op_dat dat);
void op_perf_comms(void *k_i, int nargs, op_arg *args);
#endif

void comm_timing_output();
void timing_breakdown_output(op_comm_kernel *kernel_entry);
