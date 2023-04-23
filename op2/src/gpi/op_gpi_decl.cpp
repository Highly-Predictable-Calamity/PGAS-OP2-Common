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



/* ------------------------------------------------------- */
/*      BACKEND MEMORY FUNCTIONS FOR GPI SEGMENT HEAP      */
/*         uses heap implementation taken from:            */
/*            github.com/NickWalker1/PunchOS               */
/* ------------------------------------------------------- */


MemorySegmentHeader_t *eeh_seg_start;
MemorySegmentHeader_t *enh_seg_start;
MemorySegmentHeader_t *ieh_seg_start;
MemorySegmentHeader_t *inh_seg_start;

MemorySegmentHeader_t *intialise_gpi_heap_segment(gaspi_segment_id_t seg_id, int size);

/* Initialises all segments for heap gpi dat data.*/
void op_gpi_setup_segments_heap(){
    intialise_gpi_heap_segment(EEH_HEAP_SEGMENT_ID, GPI_HEAP_SIZE);
    intialise_gpi_heap_segment(ENH_HEAP_SEGMENT_ID, GPI_HEAP_SIZE);
    intialise_gpi_heap_segment(IEH_HEAP_SEGMENT_ID, GPI_HEAP_SIZE);
    intialise_gpi_heap_segment(INH_HEAP_SEGMENT_ID, GPI_HEAP_SIZE);
}


/* Should not be called by the programmer. Use malloc or shr_malloc instead */
gaspi_offset_t op_gpi_segment_malloc(gaspi_segment_id_t  seg_id,int size){
    MemorySegmentHeader_t *currSeg;
    char *base_addr;
    /* Lookup start segment for each heap*/
    switch(seg_id){
        case EEH_HEAP_SEGMENT_ID:
            currSeg=eeh_seg_start;
            break;
        case ENH_HEAP_SEGMENT_ID:
            currSeg=enh_seg_start;
            break;
        case IEH_HEAP_SEGMENT_ID:
            currSeg=ieh_seg_start;
            break;
        case INH_HEAP_SEGMENT_ID:
            currSeg=inh_seg_start;
            break;
        default:
            GPI_FAIL("Invalid segment ID for heap allocation.\n")
    }

    base_addr=(char*)(currSeg+1);

    //traverse linked list to find one that meets conditions
    //if at end return NULL
    while(!currSeg->free || size > currSeg->size){
        currSeg=currSeg->next;
        if(!currSeg) return (gaspi_offset_t)-1;
    }
    /* Update the segment info */
    uint32_t init_size =currSeg->size;
    currSeg->size=size;
    currSeg->free=false;

    /* Quick check for heap corruption, as this should never exist that two free states are next to each other */
    if(currSeg->next && currSeg->next->free)GPI_FAIL("Heap corruption. (Likely erronous writes)");


    /* Check if we can fit a new header in + some bytes to ensure it's actually useful. */
    if(init_size>sizeof(MemorySegmentHeader_t)+64){
        MemorySegmentHeader_t* newSegment = currSeg+1;/* this should work with c struct pointer arithmetic*/
        newSegment->free=true;
        newSegment->magic=segment_magic;
        newSegment->previous=currSeg;
        newSegment->next=currSeg->next;
        newSegment->size=init_size-sizeof(MemorySegmentHeader_t)-size;

        if(currSeg->next) currSeg->next->previous=newSegment;
        currSeg->next=newSegment;

        //Make sure to return the start of the free memory space, not the space containing
        //the header information.
        return (gaspi_offset_t)(((char*)(++currSeg)) - base_addr);
    }


    //cannot fit a new segmentHeader in.

    //adjust currSeg to be initial size
    currSeg->size=init_size;
    
    //Make sure to return the start of the free memory space, not the space containing
    //the header information.
    return (gaspi_offset_t)(((char*)(++currSeg)) - base_addr);
}

/* Free the associated memory segment with addr */
void segment_free(gaspi_segment_id_t seg_id,gaspi_offset_t offset){
    
    char *base_addr;
    /* Obtain base address from segment ID*/
    switch(seg_id){
        case EEH_HEAP_SEGMENT_ID:
            base_addr=eeh_segment_ptr;
            break;
        case ENH_HEAP_SEGMENT_ID:
            base_addr=enh_segment_ptr;
            break;
        case IEH_HEAP_SEGMENT_ID:
            base_addr=ieh_segment_ptr;
            break;
        case INH_HEAP_SEGMENT_ID:
            base_addr=inh_segment_ptr;
            break;
        default: 
            GPI_FAIL("Invalid segment ID to free from.\n");
    }
    char *addr=base_addr+(unsigned int)offset; /* lossy cast but likely okay as segments likely <4GiB*/

    //Assumption made that addr is base of the free space
    MemorySegmentHeader_t *currSeg = (MemorySegmentHeader_t*) (addr - sizeof(MemorySegmentHeader_t));
    
    if(currSeg->magic!=segment_magic){
        GPI_FAIL("Attempted to free non-base address");
        return;
    }

    //sanity check
    if(currSeg->free) return;


    currSeg->free=true;

    //There can be at most one free block on either side. 
    //Eg: |Free|this block|Free|another block...
    // the case |Free|this block|Free|Free|... should never exist 

    //if the previous segment exists and is free
    if(currSeg->previous && currSeg->previous->free){
        MemorySegmentHeader_t* prevSeg = currSeg->previous;

        prevSeg->size=prevSeg->size + currSeg->size + sizeof(MemorySegmentHeader_t);        


        //update pointers to remove currSeg from linked list.
        prevSeg->next=currSeg->next;
        if(currSeg->next)
            currSeg->next->previous=prevSeg;

        //for the second if
        currSeg=prevSeg;
    }

    //if the next segment exists and is free
    if(currSeg->next && currSeg->next->free){
        currSeg->size=currSeg->size+currSeg->next->size + sizeof(MemorySegmentHeader_t);

        //update pointers to remove the next segment

        //if there is a following one
        if(currSeg->next->next){
            currSeg->next=currSeg->next->next;
            currSeg->next->previous=currSeg;
        }
        else{
            currSeg->next=0;
        }
        
    }
}


/* Allocates the segment memory, binds it to GPI segment, registers with other processes, initialises heap*/
MemorySegmentHeader_t *intialise_gpi_heap_segment(gaspi_segment_id_t seg_id, int size){
    MemorySegmentHeader_t **seg;
    char **seg_ptr;
    
    switch(seg_id){
        case EEH_HEAP_SEGMENT_ID:
            seg=&eeh_seg_start;
            seg_ptr=&eeh_heap_segment_ptr;
            break;
        case ENH_HEAP_SEGMENT_ID:
            seg=&enh_seg_start;
            seg_ptr=&enh_heap_segment_ptr;
            break;
        case IEH_HEAP_SEGMENT_ID:
            seg=&ieh_seg_start;
            seg_ptr=&ieh_heap_segment_ptr;
            break;
        case INH_HEAP_SEGMENT_ID:
            seg=&inh_seg_start;
            seg_ptr=&inh_heap_segment_ptr;
            break;
        default: 
            GPI_FAIL("Invalid segment ID to initialise.\n")
    }

    /* allocate the segment */
    *seg_ptr = (char*) xmalloc(size);
    
    int ret = gaspi_segment_use(seg_id, (gaspi_pointer_t)*seg_ptr, size, OP_GPI_WORLD, GPI_TIMEOUT, GASPI_ALLOC_DEFAULT);
    if(ret!=0){
        GPI_FAIL("Unable to use segment.\n")
    }

    *seg = (MemorySegmentHeader_t*)seg_ptr;
    (*seg)->size=size;
    (*seg)->magic=segment_magic;
    (*seg)->next=NULL;
    (*seg)->previous=NULL;
    (*seg)->free=true;


    return NULL;
}


gaspi_offset_t segment_malloc(gaspi_segment_id_t seg, int bytes){

}