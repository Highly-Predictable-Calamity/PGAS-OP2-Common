// not entirely sure if anything really needs to go in here?
//#include <op_gpi_core.h>
#include <op_perf_common.h>
#include <op_lib_core.h>
#include <mpi.h>

void op_gpi_timing_output_core() {
  if (OP_kern_max > 0) {
    if (op_is_root())
      printf("\n  count   plan time     GPI time(std)        time(std)         "
             "  GB/s      GB/s   kernel name ");
    if (op_is_root())
      printf("\n "
             "-----------------------------------------------------------------"
             "--------------------------\n");
    for (int n = 0; n < OP_kern_max; n++) {
      if (OP_kernels[n].count > 0) {
        if (OP_kernels[n].ntimes == 1 && OP_kernels[n].times[0] == 0.0f &&
            OP_kernels[n].time != 0.0f) {
          // This library is being used by an OP2 translation made with the
          // older
          // translator with older timing logic. Adjust to new logic:
          OP_kernels[n].times[0] = OP_kernels[n].time;
        }

        double kern_time = OP_kernels[n].times[0];
        for (int i=1; i<OP_kernels[n].ntimes; i++) {
          if (OP_kernels[n].times[i] > kern_time)
            kern_time = OP_kernels[n].times[i];
        }

        double moments_mpi_time[2];
        double moments_time[2];
        double moments_gpi_time[2];
        op_compute_moment_across_times(OP_kernels[n].times, OP_kernels[n].ntimes, true, &moments_time[0],
                          &moments_time[1]);
        op_compute_moment(OP_kernels[n].mpi_time, &moments_mpi_time[0],
                          &moments_mpi_time[1]);
        op_compute_moment(OP_kernels[n].gpi_time, &moments_gpi_time[0],
                          &moments_gpi_time[1]);
        if (OP_kernels[n].transfer2 < 1e-8f) {
          float transfer =
              MAX(0.0f, OP_kernels[n].transfer / (1e9f * kern_time -
                                                  OP_kernels[n].gpi_time));

          if (op_is_root())
            printf(" %6d;  %8.4f;  %8.4f(%8.4f);  %8.4f(%8.4f);  %8.4f;        "
                   " ;   %s \n",
                   OP_kernels[n].count, OP_kernels[n].plan_time,
                   moments_gpi_time[0],
                   sqrt(moments_gpi_time[1] -
                        moments_gpi_time[0] * moments_gpi_time[0]),
                   moments_time[0],
                   sqrt(moments_time[1] - moments_time[0] * moments_time[0]),
                   transfer, OP_kernels[n].name);
        } else {
          float transfer =
              MAX(0.0f, OP_kernels[n].transfer /
                            (1e9f * kern_time -
                             OP_kernels[n].plan_time - OP_kernels[n].gpi_time));
          float transfer2 =
              MAX(0.0f, OP_kernels[n].transfer2 /
                            (1e9f * kern_time -
                             OP_kernels[n].plan_time - OP_kernels[n].gpi_time));
          if (op_is_root())
            printf(" %6d;  %8.4f;  %8.4f(%8.4f);  %8.4f(%8.4f); %8.4f; %8.4f;  "
                   " %s \n",
                   OP_kernels[n].count, OP_kernels[n].plan_time,
                   moments_mpi_time[0],
                   sqrt(moments_gpi_time[1] -
                        moments_gpi_time[0] * moments_gpi_time[0]),
                   moments_time[0],
                   sqrt(moments_time[1] - moments_time[0] * moments_time[0]),
                   transfer, transfer2, OP_kernels[n].name);
        }
      }
    }
  }
}


void gpi_timing_output() {
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



    op_comm_kernel *k;
		int i = 0;
		// Print the rank-local performance information within a lockstep for each rank
		while (i < comm_size) {
			MPI_Barrier(OP_MPI_IO_WORLD);
			if (i == my_rank) {

				printf("g__________________________________________________\n");
				printf("Performance information on rank %d\n", my_rank);
				printf("Kernel        Count  total time(sec)  Avg time(sec)  \n");
				for (k = op_comm_kernel_tab; k != NULL; k = (op_comm_kernel *)k->hh.next) {
					if (k->count > 0) {
						printf("%-10s  %6d       %10.4f      %10.4f    \n", k->name, k->count,
									 k->time, k->time / k->count);
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
			} /* if i==my_rank */
			i++;
			MPI_Barrier(OP_MPI_IO_WORLD);
		}

		// Wait for all ranks to print before printing the general summary
		MPI_Barrier(OP_MPI_IO_WORLD);
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

void op_gpi_timing_output() {
  printf("op_gpi_timing_output\n");
  double max_plan_time = 0.0;
  MPI_Reduce(&OP_plan_time, &max_plan_time, 1, MPI_DOUBLE, MPI_MAX, 0,
             OP_MPI_WORLD);
  op_gpi_timing_output_core();
  if (op_is_root())
    printf("Total plan time: %8.4f\n", OP_plan_time);

  comm_timing_output();
}

