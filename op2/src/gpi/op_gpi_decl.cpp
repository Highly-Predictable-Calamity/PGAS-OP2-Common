#include <op_lib_gpi.h>
#include <op_gpi_core.h>
#include <op_util.h>
#include "gpi_utils.h"

gaspi_size_t eeh_size, enh_size, ieh_size, inh_size;

extern halo_list *OP_export_exec_list; // EEH list
extern halo_list *OP_import_exec_list; // IEH list

extern halo_list *OP_import_nonexec_list; // INH list
extern halo_list *OP_export_nonexec_list; // ENH list



int op_gpi_buffer_setup(op_dat dat, int flags){
    //Check flags
    if(!(flags & GPI_STD_DAT || flags & GPI_HEAP_DAT)){
        GPI_FAIL("Must initialise dat GPI buffer on either heap or static\n");
    }


    /* Grab the lists */
    halo_list imp_exec_list = OP_import_exec_list[dat->set->index];
    halo_list imp_nonexec_list = OP_import_nonexec_list[dat->set->index];
    
    halo_list exp_exec_list = OP_export_exec_list[dat->set->index];
    halo_list exp_nonexec_list = OP_export_nonexec_list[dat->set->index];


    /* Allocate gpi buffer */
    op_gpi_buffer gpi_buf = (op_gpi_buffer)xmalloc(sizeof(op_gpi_buffer_core));
    
    /* allocate the remote segment offsets arrays for the import lists */
    gpi_buf->remote_exec_offsets = (gaspi_offset_t*)xmalloc(sizeof(gaspi_offset_t)*exp_exec_list->ranks_size);
    gpi_buf->remote_nonexec_offsets = (gaspi_offset_t*)xmalloc(sizeof(gaspi_offset_t)*exp_nonexec_list->ranks_size);


   
    /* used to calculate offset for each rank inside the dat.
    * can use the size calculation to give the end of the previous dat.
    */
    gaspi_offset_t exec_dat_rank_offset;
    gaspi_offset_t nonexec_dat_rank_offset;

    /* Update the dat to state where the dat data starts within the segment
    *  this is used for the sending process.
    */
    if(flags & GPI_STD_DAT){
        dat->loc_eeh_seg_off=(int)eeh_size;
        dat->loc_enh_seg_off=(int)enh_size;

        exec_dat_rank_offset = ieh_size;
        nonexec_dat_rank_offset = inh_size;
    
    }
    else{ /* GPI_HEAP_DAT */

    }



    /* You receive from the import halos */
    gpi_buf->exec_recv_count = imp_exec_list->ranks_size;
    gpi_buf->nonexec_recv_count = imp_nonexec_list->ranks_size;


    /* Allocate memory for recv objs */
    gpi_buf->exec_recv_objs = (op_gpi_recv_obj *)xmalloc(gpi_buf->exec_recv_count * sizeof(op_gpi_recv_obj));
    gpi_buf->nonexec_recv_objs = (op_gpi_recv_obj *)xmalloc(gpi_buf->nonexec_recv_count * sizeof(op_gpi_recv_obj));

    dat->gpi_buffer = (void *)gpi_buf;

    /* populate the recv obj info */

    op_gpi_recv_obj *recv_obj;

    /* base offset into the exec segment of the data array */
    int exec_init = dat->set->size * dat->size;


    /* populate the exec recv objects info*/
    for (int i = 0; i < gpi_buf->exec_recv_count; i++)
    {
    recv_obj = &gpi_buf->exec_recv_objs[i];

    recv_obj->remote_rank = imp_exec_list->ranks[i];
    recv_obj->size = dat->size * imp_exec_list->sizes[i];

    /* Address calc reused from op_mpi_rt_support: op_exchange_halo */
    recv_obj->memcpy_addr = &dat->data[exec_init + imp_exec_list->disps[i] * dat->size];

    recv_obj->segment_recv_offset = exec_dat_rank_offset;

    // increment the segment offset by the size of the recv object data
    exec_dat_rank_offset += recv_obj->size;
    }

    int nonexec_init = (dat->set->size + imp_exec_list->size) * dat->size;

    /* same as above but for non-execute elements*/
    for (int i = 0; i < gpi_buf->nonexec_recv_count; i++)
    {
    recv_obj = &gpi_buf->nonexec_recv_objs[i];

    recv_obj->remote_rank = imp_nonexec_list->ranks[i];
    recv_obj->size = dat->size * imp_nonexec_list->sizes[i];

    /* Address calc reused from op_mpi_rt_support: op_exchange_halo */
    recv_obj->memcpy_addr = &dat->data[nonexec_init + imp_nonexec_list->disps[i] * dat->size];

    recv_obj->segment_recv_offset = nonexec_dat_rank_offset;

    // increment the segment offset by the size of the recv object data
    nonexec_dat_rank_offset += recv_obj->size;
    }


    // Increment by the number of elements in the halo multiplied by 
    /* Increment the reqiured segment sizes count.
    * Do not move earlier, as previous value needed for segment offset calculations.
    */
    if(flags & GPI_STD_DAT){
        eeh_size += (gaspi_size_t)exp_exec_list->size * dat->size;
        enh_size += (gaspi_size_t)exp_nonexec_list->size * dat->size;
        ieh_size += (gaspi_size_t)imp_exec_list->size * dat->size;
        inh_size += (gaspi_size_t)imp_nonexec_list->size * dat->size;
    }

    /* Populate the remote_segment_offset array in the halo_list struct.
    * This is to know where to send data to by communicating with neighbours...
    */

    /* We now know where we'd like to receieve everything, so just need to tell everyone else that...
    * Can do this very inefficiently with a LOT of sends/receives
    * I.e. one for each rank for each dat...
    * Importantly, operations with same O() communication complexity are also done in previous steps...
    */

    /* SEND */

    // Firstly need to allocate enough room for the sending MPI_Requests...
    gpi_buf->pre_exchange_hndl_s = (MPI_Request *)xmalloc(sizeof(MPI_Request) * (imp_exec_list->ranks_size + imp_nonexec_list->ranks_size));

    bool send_okay=true;

    /* Execute elements */
    for (int i = 0; i < imp_exec_list->ranks_size; i++)
    {
    recv_obj = &gpi_buf->exec_recv_objs[i];

    send_okay = send_okay &
        MPI_Isend(
        &recv_obj->segment_recv_offset,
        1,
        MPI_UNSIGNED_LONG,
        recv_obj->remote_rank,
        dat->index,
        OP_MPI_WORLD,
        &gpi_buf->pre_exchange_hndl_s[i])==MPI_SUCCESS;
    }

    /* Non-execute elements */
    for (int i = 0; i < imp_nonexec_list->ranks_size; i++)
    {
    recv_obj = &gpi_buf->nonexec_recv_objs[i];

    send_okay = send_okay &
        MPI_Isend(
        &recv_obj->segment_recv_offset,
        1,
        MPI_UNSIGNED_LONG,
        recv_obj->remote_rank,
        1 << 20 | dat->index, /* 1 in MSB -1 to indicate non-execute */
        OP_MPI_WORLD,
        &gpi_buf->pre_exchange_hndl_s[i + gpi_buf->exec_recv_count] // as sharing a MPI_Request array
        )==MPI_SUCCESS;
    }
    if(!send_okay){
    GPI_FAIL("Status code on MPI_IRecv non-zero\n");
    }

    /* RECEIVE */

    /* Receiving involves telling MPI where to put the data it recieves.
    * Data can be uniquely identified by the sending host and the tag.
    * As part of the tag, a 1 is added to the second most significant bit to indicate
    * this is for the non-execute elements */

    /* Execute elements */

    /* RECEIVE - BLOCKING */
    bool recv_okay = true;
    for (int i = 0; i < exp_exec_list->ranks_size; i++)
    {
    recv_okay = recv_okay &
        MPI_Recv(&(gpi_buf->remote_exec_offsets[i]),
                1,
                MPI_UNSIGNED_LONG,
                exp_exec_list->ranks[i],
                dat->index,
                OP_MPI_WORLD,
                MPI_STATUS_IGNORE)==MPI_SUCCESS;
    }
    for (int i = 0; i < exp_nonexec_list->ranks_size; i++)
    {
    recv_okay = recv_okay &
        MPI_Recv(&(gpi_buf->remote_nonexec_offsets[i]),
                1,
                MPI_UNSIGNED_LONG,
                exp_nonexec_list->ranks[i],
                1 << 20 | dat->index, /* 1 in MSB -1 to indicate non-exec */
                OP_MPI_WORLD,
                MPI_STATUS_IGNORE)==MPI_SUCCESS;
    }
    if(!recv_okay){
    GPI_FAIL("Status code on MPI_IRecv non-zero\n");
    }

    printf("Setup GPI stuff for %s, with buff %p\n",dat->name, dat->gpi_buffer);
    fflush(stdout);
    return 0;
}
