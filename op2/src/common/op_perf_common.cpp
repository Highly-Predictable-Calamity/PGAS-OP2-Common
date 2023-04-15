#include <op_perf_common.h>
#include <op_util.h>

/* table holding communication performance of each loop
   (accessed via a hash of loop name) */
op_comm_kernel *op_comm_kernel_tab = NULL;

op_comm_kernel *get_kernel_entry(const char *name){
  op_comm_kernel *kernel_entry;

  HASH_FIND_STR(op_comm_kernel_tab, name, kernel_entry);
  if (kernel_entry == NULL) {
    kernel_entry = (op_comm_kernel *)xmalloc(sizeof(op_comm_kernel));
    kernel_entry->num_indices = 0;
    kernel_entry->time = 0.0;
    kernel_entry->count = 0;
    kernel_entry->breakdown_cap = 0;
    kernel_entry->breakdown = NULL;
    strncpy((char *)kernel_entry->name, name, NAMESIZE);
    HASH_ADD_STR(op_comm_kernel_tab, name, kernel_entry);
  }

  return kernel_entry;
}

/*******************************************************************************
 * Routine to measure timing for an op_par_loop / kernel
 *******************************************************************************/
void *op_comm_perf_time(const char *name, double time) {
  op_comm_kernel *kernel_entry;
  kernel_entry = get_kernel_entry(name);
  kernel_entry->count += 1;
  kernel_entry->time += time;
  return (void *)kernel_entry;
}

kernel_breakdown *get_breakdown(op_comm_kernel *kernel_entry, const char *breakdown_name){
  for (int i=0; i<kernel_entry->breakdown_cap; i++){
    if (strcmp(breakdown_name, kernel_entry->breakdown[i].name) == 0){
      return &kernel_entry->breakdown[i];
    }
  }
  kernel_entry->breakdown_cap += 1;
  kernel_entry->breakdown = (kernel_breakdown *)xrealloc(kernel_entry->breakdown, sizeof(kernel_breakdown)*kernel_entry->breakdown_cap);
  
  kernel_breakdown *new_breakdown = kernel_entry->breakdown + kernel_entry->breakdown_cap-1;
  new_breakdown->count = 0;
  new_breakdown->time = 0.0;
  strncpy((char*)new_breakdown->name, breakdown_name, NAMESIZE);
  return new_breakdown;
}

void *op_comm_perf_time_breakdown(const char *kernel_name, const char *breakdown_name, double time){
  op_comm_kernel *kernel_entry;
  kernel_entry = get_kernel_entry(kernel_name);
  kernel_breakdown *breakdown = get_breakdown(kernel_entry, breakdown_name);
  breakdown->count += 1;
  breakdown->time += time;
  return NULL;
}


#ifdef COMM_PERF // not entirely sure what this bit is for, but mpi references have been removed

/*******************************************************************************
 * Routine to linear search comm_info array in an op_comm_kernel for an op_dat
 *******************************************************************************/
int search_op_comm_kernel(op_dat dat, op_comm_kernel *kernal, int num_indices) {
  for (int i = 0; i < num_indices; i++) {
    if (strcmp((kernal->comm_info[i])->name, dat->name) == 0 &&
        (kernal->comm_info[i])->size == dat->size) {
      return i;
    }
  }

  return -1;
}

/*******************************************************************************
 * Routine to measure message sizes exchanged in an op_par_loop / kernel
 *******************************************************************************/
void op__perf_comm(void *k_i, op_dat dat) {
  halo_list exp_exec_list = OP_export_exec_list[dat->set->index];
  halo_list exp_nonexec_list = OP_export_nonexec_list[dat->set->index];
  int tot_halo_size =
      (exp_exec_list->size + exp_nonexec_list->size) * (size_t)dat->size;

  op_comm_kernel *kernel_entry = (op_comm_kernel *)k_i;
  int num_indices = kernel_entry->num_indices;

  if (num_indices == 0) {
    // set capcity of comm_info array
    kernel_entry->cap = 20;
    op_dat_comm_info dat_comm =
        (op_dat_comm_info)xmalloc(sizeof(op_dat_comm_info_core));
    kernel_entry->comm_info = (op_dat_comm_info *)xmalloc(
        sizeof(op_dat_comm_info *) * (kernel_entry->cap));
    strncpy((char *)dat_comm->name, dat->name, 20);
    dat_comm->size = dat->size;
    dat_comm->index = dat->index;
    dat_comm->count = 0;
    dat_comm->bytes = 0;

    // add first values
    dat_comm->count += 1;
    dat_comm->bytes += tot_halo_size;

    kernel_entry->comm_info[num_indices] = dat_comm;
    kernel_entry->num_indices++;
  } else {
    int index = search_op_comm_kernel(dat, kernel_entry, num_indices);
    if (index < 0) {
      // increase capacity of comm_info array
      if (num_indices >= kernel_entry->cap) {
        kernel_entry->cap = kernel_entry->cap * 2;
        kernel_entry->comm_info = (op_dat_comm_info *)xrealloc(
            kernel_entry->comm_info,
            sizeof(op_dat_comm_info *) * (kernel_entry->cap));
      }

      op_dat_comm_info dat_comm =
          (op_dat_comm_info)xmalloc(sizeof(op_dat_comm_info_core));

      strncpy((char *)dat_comm->name, dat->name, 20);
      dat_comm->size = dat->size;
      dat_comm->index = dat->index;
      dat_comm->count = 0;
      dat_comm->bytes = 0;

      // add first values
      dat_comm->count += 1;
      dat_comm->bytes += tot_halo_size;

      kernel_entry->comm_info[num_indices] = dat_comm;
      kernel_entry->num_indices++;
    } else {
      kernel_entry->comm_info[index]->count += 1;
      kernel_entry->comm_info[index]->bytes += tot_halo_size;
    }
  }
}
#endif

#ifdef COMM_PERF
void op_perf_comms(void *k_i, int nargs, op_arg *args) {

  for (int n = 0; n < nargs; n++) {
    if (args[n].argtype == OP_ARG_DAT && args[n].sent == 2) {
      op_perf_comm(k_i, (&args[n])->dat);
    }
  }
}
#endif

void comm_timing_output() {
  int my_rank, comm_size;
  MPI_Comm OP_MPI_IO_WORLD;
  MPI_Comm_dup(OP_MPI_WORLD, &OP_MPI_IO_WORLD);
  MPI_Comm_rank(OP_MPI_IO_WORLD, &my_rank);
  MPI_Comm_size(OP_MPI_IO_WORLD, &comm_size);

  unsigned int count, tot_count;
  count = HASH_COUNT(op_comm_kernel_tab);
  MPI_Allreduce(&count, &tot_count, 1, MPI_INT, MPI_SUM, OP_MPI_IO_WORLD);

  if (tot_count > 0) {
    double tot_time;
    double avg_time;

    printf("___________________________________________________\n");
    printf("Performance information on rank %d\n", my_rank);
    printf("Kernel        Count  total time(sec)  Avg time(sec)  \n");

    op_comm_kernel *k;
    for (k = op_comm_kernel_tab; k != NULL; k = (op_comm_kernel *)k->hh.next) {
      if (k->count > 0) {
        printf("%-10s  %6d       %10.4f      %10.4f    \n", k->name, k->count,
               k->time, k->time / k->count);
        fflush(stdout);
        timing_breakdown_output(k);

#ifdef COMM_PERF
        if (k->num_indices > 0) {
          printf("halo exchanges:  ");
          for (int i = 0; i < k->num_indices; i++)
            printf("%10s ", k->comm_info[i]->name);
          printf("\n");
          printf("       count  :  ");
          for (int i = 0; i < k->num_indices; i++)
            printf("%10d ", k->comm_info[i]->count);
          printf("\n");
          printf("total(Kbytes) :  ");
          for (int i = 0; i < k->num_indices; i++)
            printf("%10d ", k->comm_info[i]->bytes / 1024);
          printf("\n");
          printf("average(bytes):  ");
          for (int i = 0; i < k->num_indices; i++)
            printf("%10d ", k->comm_info[i]->bytes / k->comm_info[i]->count);
          printf("\n");
        } else {
          printf("halo exchanges:  %10s\n", "NONE");
        }
        printf("---------------------------------------------------\n");
#endif
      }
    }
    printf("___________________________________________________\n");

    if (my_rank == MPI_ROOT) {
      printf("___________________________________________________\n");
      printf("\nKernel        Count   Max time(sec)   Avg time(sec)  \n");
    }

    for (k = op_comm_kernel_tab; k != NULL; k = (op_comm_kernel *)k->hh.next) {
      MPI_Reduce(&(k->count), &count, 1, MPI_INT, MPI_MAX, MPI_ROOT,
                 OP_MPI_IO_WORLD);
      MPI_Reduce(&(k->time), &avg_time, 1, MPI_DOUBLE, MPI_SUM, MPI_ROOT,
                 OP_MPI_IO_WORLD);
      MPI_Reduce(&(k->time), &tot_time, 1, MPI_DOUBLE, MPI_MAX, MPI_ROOT,
                 OP_MPI_IO_WORLD);

      if (my_rank == MPI_ROOT && count > 0) {
        printf("%-10s  %6d       %10.4f      %10.4f    \n", k->name, count,
               tot_time, (avg_time) / comm_size);
      }
      tot_time = avg_time = 0.0;
    }
  }
  MPI_Comm_free(&OP_MPI_IO_WORLD);
}

void timing_breakdown_output(op_comm_kernel *kernel_entry){
  for (int i=0; i<kernel_entry->breakdown_cap; i++){
    kernel_breakdown bd = kernel_entry->breakdown[i];
    const char* line = (i==kernel_entry->breakdown_cap-1) ? "└─" : "├─";
    printf("%s ", line);
    printf("%-7s  %6d       %10.4f      %10.4f    \n",
           bd.name, bd.count, bd.time, bd.time / bd.count);
    if (i==kernel_entry->breakdown_cap-1) printf("\n");
  }
}