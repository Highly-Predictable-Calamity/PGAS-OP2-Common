#include <GASPI.h>

#include "op_mpi_core.h" //Include the mpi stuff as it's the bigger thing


#include "../src/gpi/gpi_utils.h"

#define EEH_SEGMENT_ID 1
#define ENH_SEGMENT_ID 2
#define IEH_SEGMENT_ID 3
#define INH_SEGMENT_ID 4

#define MSC_SEGMENT_ID 5

extern gaspi_group_t OP_GPI_WORLD;
extern gaspi_group_t OP_GPI_GLOBAL;

#define OP2_GPI_QUEUE_ID 1 /* All remote comm. uses the same queue ID for simplicity. */

/* Given in chars as offsets stored as bytes so pointer arithmetic correct and easy*/
extern char *eeh_segment_ptr;
extern char *ieh_segment_ptr;
extern char *enh_segment_ptr;
extern char *inh_segment_ptr;

extern char *msc_segment_ptr;

/* Struct storing information regarding the expected dat elements from who, where, and where to copy to */
typedef struct{
  gaspi_rank_t        remote_rank; /* Rank receiving from - (used to identify the struct) */
  gaspi_offset_t      segment_recv_offset; /* segment offset in bytes for the receieved information (where the remote write landed) */
  void*               memcpy_addr; /* Where to memcpy the received data to. I.e. the vaddr of the location inside the data array */
  int                 size; /* Number of bytes */
/*?smart linked list entry struct?*/
} op_gpi_recv_obj; 

struct op_gpi_buffer_core{
  int exec_recv_count; /* Number of recieves for import execute segment expect (i.e. number of remote ranks)*/
  int nonexec_recv_count; /* Number of recieves for import non-execute segment expect (i.e number of remote ranks)*/
  op_gpi_recv_obj *exec_recv_objs; /*  For exec elements of this dat, one for each of the expected notifications*/
  op_gpi_recv_obj *nonexec_recv_objs; /* For nonexec elements of this dat , one for each of the expected notifications*/
  MPI_Request *pre_exchange_hndl_s; /* data pre exchange handles for sends */
  MPI_Request *pre_exchange_hndl_r; /* data pre exchange handles for receives */
};

typedef op_gpi_buffer_core *op_gpi_buffer;

void op_gpi_exchange_halo(op_arg *arg, int exec_flag);
void op_gpi_exchange_halo_partial(op_arg *arg, int exec_flag);

void op_gpi_waitall(op_arg *arg);
void op_gpi_waitall_args(int nargs, op_arg *args);