#include <GASPI.h>

#include <op_lib_c.h>
#include <op_lib_core.h>
#include <op_util.h>

#include <op_lib_gpi.h>
#include <op_lib_mpi.h>
#include <op_mpi_core.h>

#include <op_rt_support.h> /* Check this- likely not needed.*/

#include <op_gpi_core.h> 
#include <op_perf_common.h>

#include "gpi_utils.h"


/* GPI reimplementation of op_exchange_halo originally found in op_mpi_rt_support.cpp 
 * IS_COMMON 
 * Lots of this is common, so can be put there. 
*/
void op_gpi_exchange_halo(op_arg *arg, int exec_flag){
    op_dat dat = arg->dat;
    op_gpi_buffer gpi_buf = (op_gpi_buffer) dat->gpi_buffer;

    //If it's not in use, don't bother!
    if(arg->opt ==0)
        return;

    //Check if arg already sent
    if(arg->sent ==1){
        GPI_FAIL("Error: halo exchange already in flight for dat %s\n",dat->name);
    }

    // For a directly accessed op_dat do not do halo exchanges if not executing
    // over
    // redundant compute block
    if(exec_flag == 0 && arg->idx == -1)
        return;

    arg->sent =0; //reset flag (TODO seems unneccesary but keep anyway)

    // need to exchange both direct and indirect data sets if they're dirty
    // return if not R/RW or if not dirty.
    if(!(arg->acc == OP_READ || arg->acc == OP_RW) 
        || (dat->dirtybit !=1))
            return;
    

    //Grab the halo lists
    halo_list imp_exec_list = OP_import_exec_list[dat->set->index];
    halo_list imp_nonexec_list = OP_import_nonexec_list[dat->set->index];

    halo_list exp_exec_list = OP_export_exec_list[dat->set->index];
    halo_list exp_nonexec_list = OP_export_nonexec_list[dat->set->index];

    int gpi_rank;
    gaspi_proc_rank((gaspi_rank_t*)&gpi_rank);

    //-------first exchange exec elements related to this data array--------

    //sanity checks
    if (compare_sets(imp_exec_list->set, dat->set) == 0 ){
        GPI_FAIL("Import list and set mismatch\n");
    }
    if (compare_sets(exp_exec_list->set, dat->set) == 0 ){
        GPI_FAIL("Export list and set mismatch\n");
    }

    //dat offset inside eeh segment
    //Note - changed to int to not perform addition on pointer type
    // and simplifies offset logic as operations are now performed on byte count.
    void *dat_offset_addr; 
    if(gpi_buf->is_dynamic){
        dat_offset_addr =(void*)(eeh_heap_segment_ptr + dat->loc_eeh_seg_off);
    }
    else{
        dat_offset_addr = (void*)(eeh_segment_ptr + dat->loc_eeh_seg_off);
    }

    int set_elem_index;
    for (int i = 0; i < exp_exec_list->ranks_size; i++) {
      for (int j = 0; j < exp_exec_list->sizes[i]; j++) {

        set_elem_index = exp_exec_list->list[exp_exec_list->disps[i] + j];
        //Can reuse the exp_exec_list->disps[i] as this gives the per rank displacement into the dat buffer.

        //memcpy into eeh segment appropriately
        memcpy(((void*)dat_offset_addr + exp_exec_list->disps[i]* dat->size + j* dat->size),
               (void *)&dat->data[dat->size * (set_elem_index)], dat->size);
      }

      //get remote offsets for that rank
      gaspi_offset_t remote_exec_offset = (gaspi_offset_t) gpi_buf->remote_exec_offsets[i];

      gaspi_offset_t local_offset = (gaspi_offset_t) dat->loc_eeh_seg_off+ exp_exec_list->disps[i]*dat->size;

        /* Wait for acknowledgement before send */
        if(gpi_buf->exec_sent_acks[i]){

            gaspi_notification_id_t wait_id;
            gaspi_notification_t wait_val;
            //printf("Waiting for acknowledgement\n");

            GPI_SAFE(gaspi_notify_waitsome(
                EEH_SEGMENT_ID + gpi_buf->is_dynamic*DYNAMIC_SEG_ID_OFFSET, /* local segment */
                dat->index<<NOTIF_SHIFT | exp_exec_list->ranks[i],
                1,
                &wait_id,
                GPI_TIMEOUT
            ) )
            
            int ack_rank = wait_id & ((1<<NOTIF_SHIFT)-1); /* Filter out the upper half */
            int ack_dat_index= (int) wait_id>> NOTIF_SHIFT; /* Get only the upper half*/
            
            //Sanity check
            if(ack_dat_index!=dat->index || ack_rank !=exp_exec_list->ranks[i]){
                GPI_FAIL("Accepted incorrect ack notification from unexpected dat.\n Expected: dat %d, rank %d got: dat %d, rank%d\n",dat->index, ack_dat_index, exp_exec_list->ranks[i],ack_rank);
            }

            GPI_SAFE(gaspi_notify_reset(
                EEH_SEGMENT_ID + gpi_buf->is_dynamic*DYNAMIC_SEG_ID_OFFSET, /* local segment */
                wait_id,
                &wait_val
            ))

            gpi_buf->exec_sent_acks[i]=0;
            //printf("acknowledgement obtained\n");
        }


      GPI_QUEUE_SAFE( gaspi_write_notify(
                        EEH_SEGMENT_ID + gpi_buf->is_dynamic*DYNAMIC_SEG_ID_OFFSET, /* local segment id*/
                        local_offset, /* local segment offset*/
                        exp_exec_list->ranks[i], /* remote rank*/
                        IEH_SEGMENT_ID  + gpi_buf->is_dynamic*DYNAMIC_SEG_ID_OFFSET, /* remote segment id*/
                        remote_exec_offset, /* remote offset*/
                        dat->size * exp_exec_list->sizes[i], /* send size*/
                        dat->index <<NOTIF_SHIFT | gpi_rank, /* notification id*/
                        1, /* notification value, 1 bit added for non-zero notif values */
                        OP2_GPI_QUEUE_ID, /* queue id*/
                        GPI_TIMEOUT /* timeout*/
                        ), OP2_GPI_QUEUE_ID )

#ifdef GPI_VERBOSE
      printf("Rank %d sent execute %s dat data to rank %d with not_ID %d \n",gpi_rank,dat->name, exp_exec_list->ranks[i],dat->index <<NOTIF_SHIFT | gpi_rank);
      fflush(stdout);
#endif
        gpi_buf->exec_sent_acks[i]=1;
    }


    //Second exchange for nonexec elements. 
    if (compare_sets(imp_nonexec_list->set, dat->set) == 0){
        GPI_FAIL("Error: Non-Import list and set mismatch");
    }

    if(compare_sets(exp_nonexec_list->set, dat->set) == 0){
        GPI_FAIL("Error: Non-Export list and set mismatch");
    }

    if(gpi_buf->is_dynamic){
        dat_offset_addr = (void*)(enh_heap_segment_ptr + dat->loc_enh_seg_off);
    }
    else{
        dat_offset_addr = (void*)(enh_segment_ptr + dat->loc_enh_seg_off);
    }

    for (int i =0; i < exp_nonexec_list->ranks_size; i++){
        for (int j=0;j<exp_nonexec_list->sizes[i];j++){
            set_elem_index = exp_nonexec_list->list[exp_nonexec_list->disps[i] + j];

            memcpy((void*)dat_offset_addr+exp_nonexec_list->disps[i]* dat->size + j *dat->size,
                   (void*)&dat->data[dat->size * (set_elem_index)],
                    dat->size);
        }
        
        gaspi_offset_t remote_nonexec_offset = (gaspi_offset_t) gpi_buf->remote_nonexec_offsets[i];

        /* Wait for acknowledgement before send */
        if(gpi_buf->nonexec_sent_acks[i]){

            gaspi_notification_id_t wait_id;
            gaspi_notification_t wait_val;
            //printf("Waiting for acknowledgement\n");

            GPI_SAFE(gaspi_notify_waitsome(
                ENH_SEGMENT_ID + gpi_buf->is_dynamic*DYNAMIC_SEG_ID_OFFSET, /* local segment */
                dat->index<<NOTIF_SHIFT | exp_nonexec_list->ranks[i],
                1,
                &wait_id,
                GPI_TIMEOUT
            ) )

            int ack_rank = wait_id & ((1<<NOTIF_SHIFT)-1); /* Filter out the upper half */
            int ack_dat_index= (int) wait_id>> NOTIF_SHIFT; /* Get only the upper half*/
            
            //Sanity check
            if(ack_dat_index!=dat->index || ack_rank !=exp_nonexec_list->ranks[i]){
                GPI_FAIL("Accepted incorrect ack notification from unexpected dat.\n Expected: dat %d, rank %d got: dat %d, rank%d\n",dat->index, ack_dat_index, exp_nonexec_list->ranks[i],ack_rank)
            }

            GPI_SAFE(gaspi_notify_reset(
                ENH_SEGMENT_ID + gpi_buf->is_dynamic*DYNAMIC_SEG_ID_OFFSET, /* local segment */
                wait_id,
                &wait_val
            ))

            gpi_buf->nonexec_sent_acks[i]=0;
            //printf("acknowledgement obtained\n");
        }


        GPI_QUEUE_SAFE( gaspi_write_notify(
                           ENH_SEGMENT_ID + gpi_buf->is_dynamic*DYNAMIC_SEG_ID_OFFSET, /* local segment */
                           (gaspi_offset_t) dat->loc_enh_seg_off + exp_nonexec_list->disps[i]*dat->size, /* local segment offset*/
                           exp_nonexec_list->ranks[i], /* remote rank*/
                           INH_SEGMENT_ID + gpi_buf->is_dynamic*DYNAMIC_SEG_ID_OFFSET, /* remote segment */
                           remote_nonexec_offset, /* remote segment offset*/
                           dat->size * exp_nonexec_list->sizes[i], /* data to send (in bytes)*/
                           dat->index<<NOTIF_SHIFT | gpi_rank, /* notification id*/
                           1, /* notification value. 1 added for non-zero notif values */
                           OP2_GPI_QUEUE_ID, /* queue id*/
                           GPI_TIMEOUT /* timeout */
                           ), OP2_GPI_QUEUE_ID )

        
        gpi_buf->nonexec_sent_acks[i]=1;

#ifdef GPI_VERBOSE
        printf("Rank %d sent non-execute %s dat data to rank %d with not_ID %d.\n",gpi_rank,dat->name, exp_nonexec_list->ranks[i], dat->index <<NOTIF_SHIFT | gpi_rank);
        fflush(stdout);
#endif
    }

    //Finish up
    dat->dirtybit =0;
    arg->sent=1;
}


/* Wait for a single arg
 * equivalent to op_mpi_waitall function
 * definitey NOT_COMMON
 */
void op_gpi_waitall(op_arg *arg){
    //Check skip conditions
    if(!(arg->opt && arg->argtype == OP_ARG_DAT && arg->sent ==1))
        return;
    
    gaspi_rank_t rank;
    gaspi_proc_rank(&rank);

    op_dat dat = arg->dat;

    op_gpi_buffer buff = (op_gpi_buffer)dat->gpi_buffer;

    op_gpi_recv_obj *exec_recv_objs = buff->exec_recv_objs;
    op_gpi_recv_obj *nonexec_recv_objs = buff->nonexec_recv_objs;


    /* Stores the op_dat_index */
    gaspi_notification_id_t notif_id;
    gaspi_notification_t    notif_value;

    int recv_rank, recv_dat_index;

#ifdef GPI_VERBOSE
    printf("Rank %d expects %d exec receives from ranks:\n",rank,buff->exec_recv_count);
    for(int i=0;i<buff->exec_recv_count;i++){
        printf("%d ",buff->exec_recv_objs[i].remote_rank);
    }
    printf("\n");
    fflush(stdout);
#endif


    /* TODO move somewhere else and store */
    int max_expected_rank=0;
    for(int i=0;i<buff->exec_recv_count;i++){
        if(buff->exec_recv_objs[i].remote_rank>max_expected_rank)
            max_expected_rank=buff->exec_recv_objs[i].remote_rank;
    }


    /* Receive for exec elements*/
    for(int i=0;i<buff->exec_recv_count;i++){
        GPI_SAFE( gaspi_notify_waitsome(IEH_SEGMENT_ID + buff->is_dynamic*DYNAMIC_SEG_ID_OFFSET, 
                            dat->index << NOTIF_SHIFT,
                            10, /* Should be max_expected rank*/
                            &notif_id, /* Notification id should be the dat index*/
                            GPI_TIMEOUT) )

#ifdef GPI_VERBOSE
        printf("Rank %d received exec not_ID %d.\n",rank,notif_id);
#endif

        recv_rank = notif_id & ((1<<NOTIF_SHIFT)-1); /* Filter out the upper half */
        recv_dat_index= (int) notif_id>> NOTIF_SHIFT; /* Get only the upper half*/
        
        //Sanity check
        if(recv_dat_index!=dat->index){
            GPI_FAIL("Accepted exec notification from unexpected dat.\n Expected: %d got:%d\n",dat->index, notif_id);
        }
        
        //store and reset notification value
        GPI_SAFE( gaspi_notify_reset(IEH_SEGMENT_ID + buff->is_dynamic*DYNAMIC_SEG_ID_OFFSET,
                            notif_id,
                            &notif_value) ) 

        /* Send acknowledgement back to sender */
        GPI_QUEUE_SAFE(gaspi_notify(
                EEH_SEGMENT_ID + buff->is_dynamic*DYNAMIC_SEG_ID_OFFSET, /* segment */
                recv_rank,
                recv_dat_index << NOTIF_SHIFT | rank,
                1,
                ACK_QUEUE,
                GPI_TIMEOUT
        ), ACK_QUEUE)



        //lookup recv object
        int obj_idx=0;
        while((obj_idx < buff->exec_recv_count) && (exec_recv_objs[obj_idx].remote_rank != recv_rank))
            obj_idx++;
        

        //check if it didn't find it...
        if(obj_idx >= buff->exec_recv_count){
            printf("Dumping exec recv objs:\n");
            printf("Exec recv count: %d\n",buff->exec_recv_count);
            for(int i=0;i<buff->exec_recv_count;i++){
              op_gpi_recv_obj *obj = &buff->exec_recv_objs[i];
              printf(" - remote rank: %d\n - segment_recv_offset: %d\n - memcpy addr: %p\n - size: %d\n\n",obj->remote_rank,obj->segment_recv_offset, obj->memcpy_addr, obj->size);
            }
            fflush(stdout);
            GPI_FAIL("Unable to find exec recv object from rank %d in dat %s\n",recv_rank, dat->name); 
        }


        //Use to memcpy data
        op_gpi_recv_obj *obj = &exec_recv_objs[obj_idx]; /* not neccessary but looks nicer later*/

        op_timers_core(&c1, &t1);
        
        char *segment_ptr = buff->is_dynamic ? ieh_heap_segment_ptr : ieh_segment_ptr;

        // Copy the data into the op_dat->data array
        memcpy(obj->memcpy_addr, (void*) (segment_ptr + obj->segment_recv_offset), obj->size);
        op_timers_core(&c2, &t2);
        op_comm_perf_time("memcpy",t2-t1);
        
        /* Send acknowledgement back to sender */
        GPI_QUEUE_SAFE(gaspi_notify(
                EEH_SEGMENT_ID + buff->is_dynamic*DYNAMIC_SEG_ID_OFFSET, /* segment */
                recv_rank,
                recv_dat_index << NOTIF_SHIFT | rank,
                1,
                ACK_QUEUE,
                GPI_TIMEOUT
        ), ACK_QUEUE)

#ifdef GPI_VERBOSE  
        printf("Rank %d successfully handled notification from rank %d for exec dat data %s.\n",rank, recv_rank,dat->name);
        fflush(stdout);
#endif
    }
    
#ifdef GPI_VERBOSE
    printf("Rank %d received neccessary exec elements for %s\n",rank,dat->name);
    fflush(stdout);

    printf("Rank %d expects %d non-exec receives from ranks:\n",rank,buff->nonexec_recv_count);
    for(int i=0;i<buff->nonexec_recv_count;i++){
        printf("%d ",buff->nonexec_recv_objs[i].remote_rank);
    }
    printf("\n");
    fflush(stdout);
#endif


    /* TODO move to store elsewhere */
    for(int i=0;i<buff->nonexec_recv_count;i++){
        if(buff->nonexec_recv_objs[i].remote_rank>max_expected_rank)
            max_expected_rank=buff->nonexec_recv_objs[i].remote_rank;
    }


    /* Receive for nonexec elements*/
    for(int i=0;i<buff->nonexec_recv_count;i++){
        GPI_SAFE( gaspi_notify_waitsome(INH_SEGMENT_ID + buff->is_dynamic*DYNAMIC_SEG_ID_OFFSET, 
                            dat->index<<NOTIF_SHIFT,
                            10, /* Should be max_expected_rank*/
                            &notif_id, /* Notification id should be the dat index*/
                            GPI_TIMEOUT) )


#ifdef GPI_VERBOSE
        printf("Rank %d received non_exec not_ID %d.\n",rank, notif_id);
#endif

        recv_rank = (int) notif_id & ((1<<NOTIF_SHIFT)-1); /* Filter out the upper half */
        recv_dat_index= (int) notif_id>>NOTIF_SHIFT; /* get only the upper half*/
        
        if(recv_dat_index!=dat->index){
            GPI_FAIL("Accepted nonexec notification from unexpected dat.\n Expected: %d got:%d\n",dat->index, notif_id);
        }
        
        //store and reset notification value
        GPI_SAFE( gaspi_notify_reset(INH_SEGMENT_ID + buff->is_dynamic*DYNAMIC_SEG_ID_OFFSET,
                            notif_id,
                            &notif_value) )
        
        /* Send acknowledgement back to sender */
        GPI_QUEUE_SAFE(gaspi_notify(
                ENH_SEGMENT_ID + buff->is_dynamic*DYNAMIC_SEG_ID_OFFSET, /* segment */
                recv_rank,
                recv_dat_index << NOTIF_SHIFT | rank,
                1,
                ACK_QUEUE,
                GPI_TIMEOUT
        ), ACK_QUEUE)





        //lookup recv object
        int obj_idx=0;
        while(obj_idx < buff->nonexec_recv_count && nonexec_recv_objs[obj_idx].remote_rank != recv_rank)
            obj_idx++;
        
        //check if it didn't find it...
        if(obj_idx >= buff->nonexec_recv_count){
            printf("Dumping exec recv objs:\n");
            printf("Exec recv count: %d\n",buff->nonexec_recv_count);
            for(int i=0;i<buff->nonexec_recv_count;i++){
              op_gpi_recv_obj *obj = &buff->nonexec_recv_objs[i];
              printf(" - remote rank: %d\n - segment_recv_offset: %d\n - memcpy addr: %p\n - size: %d\n\n",obj->remote_rank,obj->segment_recv_offset, obj->memcpy_addr, obj->size);
            }
            fflush(stdout);
            GPI_FAIL("Unable to find nonexec recv object from rank %d in dat %s\n",recv_rank, dat->name); 
        }


        //Use to memcpy data
        op_gpi_recv_obj *obj = &nonexec_recv_objs[obj_idx]; /* not neccessary but looks nicer later*/
        
        op_timers_core(&c1, &t1);
        char *segment_ptr = buff->is_dynamic ? inh_heap_segment_ptr : inh_segment_ptr;

        // Copy the data into the op_dat->data array
        memcpy(obj->memcpy_addr, (void*) (segment_ptr + obj->segment_recv_offset), obj->size);
        op_timers_core(&c2, &t2);
        op_comm_perf_time("memcpy",t2-t1);


        /* Send acknowledgement back to sender */
        GPI_QUEUE_SAFE(gaspi_notify(
                ENH_SEGMENT_ID + buff->is_dynamic*DYNAMIC_SEG_ID_OFFSET, /* segment */
                recv_rank,
                recv_dat_index << NOTIF_SHIFT | rank,
                1,
                ACK_QUEUE,
                GPI_TIMEOUT
        ), ACK_QUEUE)


#ifdef GPI_VERBOSE
        printf("Rank %d successfully handled notification from rank %d for nonexec dat data %s.\n",rank,recv_rank,dat->name);
        fflush(stdout);
#endif
    }

#ifdef GPI_VERBOSE
    printf("Rank %d receievd neccessary nonexec elements for %s\n",rank,dat->name);
    fflush(stdout);
#endif

    //Do partial halo stuff
    if(arg->map != OP_ID && OP_map_partial_exchange[arg->map->index]){
        GPI_FAIL("Not implemented partial exchange\n");
        /*
        halo_list imp_nonexec_list = OP_import_nonexec_permap[arg->map->index];
        int init = OP_export_nonexec_permap[arg->map->index]->size;
        char *buffer =
            &((op_mpi_buffer)(dat->mpi_buffer))->buf_nonexec[init * dat->size];
        for (int i = 0; i < imp_nonexec_list->size; i++) {
            int set_elem_index = imp_nonexec_list->list[i];
            memcpy((void *)&dat->data[dat->size * (set_elem_index)],
                &buffer[i * dat->size], dat->size);
        }
        */
    }
}

void op_gpi_exchange_halo_partial(op_arg *arg, int exec_flag){
    GPI_FAIL("Function is not implemented\n");
}