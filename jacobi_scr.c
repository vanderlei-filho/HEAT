#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <mpi.h>
#include <scr.h>
#include "jacobi.h"


static int rank = MPI_PROC_NULL;
static int iteration = 0;
static int verbose = 1;

// SCR query variables
static char *scr_prefix;
static int step;

// Calculate the (cumulative) execution time of SCR functions and print them at the end
//#define CKPT_DEBUG

// Single checkpoint for all processes or one checkpoint per process
//#define ONE_CKPT

// Use SCR_Need_checkpoint() or manual checkpoint check
//#define USE_SCR_NEED_CHECKPOINT

#ifdef CKPT_DEBUG
    static double t1;
    static double t2;
    static double t_have_restart     = 0.0;
    static double t_start_restart    = 0.0;
    static double t_route_file_ch    = 0.0;
    static double t_complete_restart = 0.0;
    static double t_need_checkpoint  = 0.0;
    static double t_start_output     = 0.0;
    static double t_route_file_out   = 0.0;
    static double t_complete_output  = 0.0;
#endif

/**
 * Extracts the final number from a string.
 * @param str Pointer to the string.
 * @return    The final number in the string.
*/
static int extract_final_number(char *str)
{
    int length = strlen(str);
    int number = 0;
    int multiplier = 1;

    for (int i = length - 1; i >= 0; i--)
    {
        if (isdigit(str[i]))
        {
            number += (str[i] - '0') * multiplier;
            multiplier *= 10;
        }
        else
        {
            break;
        }
    }

    return number;
}

/**
 * Reads a checkpoint file and stores the data in the input buffer.
 *
 * @param file   Pointer to the full path to a checkpoint file.
 * @param buf    Pointer to the pointer to the input buffer.
 * @param length Quantity of elements in the input buffer.
 * @return       1 if the checkpoint file was read successfully; 0 otherwise.
 */
static int read_ch(char *file, TYPE **buf, int length)
{
    int i, rc, valid = 1;
    size_t return_value;
    TYPE *read_buf;
    FILE *pFile;

    read_buf = (TYPE *)malloc(sizeof(TYPE) * length);
    if (NULL == read_buf)
    {
        printf("%d: Could not allocate memory for read_buf\n", rank);
        valid = 0;
    }
    else
    {
        pFile = fopen(file, "rb");
        if (NULL == pFile)
        {
            printf("%d: Could not open file %s\n", rank, file);
            valid = 0;
        }
        else
        {
            return_value = fread(read_buf, sizeof(TYPE), length, pFile);
            if (length != return_value)
            {
                printf("%d: Error reading %s\n", rank, file);
                valid = 0;
            }
            else
            {
                for (i = 0; i < length; i++)
                {
                    (*buf)[i] = read_buf[i];
                }
            }
            rc = fclose(pFile);
            if (0 != rc)
            {
                printf("%d: Error closing %s\n", rank, file);
                valid = 0;
            }
        }
        free(read_buf);
    }
    return valid;
}

/**
 * Checks if a checkpoint file exists and restarts from there.
 *
 * @param name   Pointer to the checkpoint filename.
 * @param buf    Pointer to the pointer to the input buffer.
 * @param length Quantity of elements in the input buffer.
 */
static void try_restart(char *name, TYPE **buf, int length)
{
    int scr_retval, have_restart, ckpt_iteration;
    int found_checkpoint = 0;
    int restarted = 0;
    char dset[SCR_MAX_FILENAME];
    char path[SCR_MAX_FILENAME];
    char file[SCR_MAX_FILENAME];
    do
    {
        if (verbose && 0 == rank)
        {
            printf("Checking for restart...\n");
        }

        #ifdef CKPT_DEBUG
            t1 = MPI_Wtime();
        #endif

        scr_retval = SCR_Have_restart(&have_restart, dset);

        #ifdef CKPT_DEBUG
            t2 = MPI_Wtime();
            t_have_restart += t2 - t1;
        #endif

        if (SCR_SUCCESS != scr_retval)
        {
            printf("%d: failed calling SCR_Have_restart: %d: @%s:%d\n",
                   rank, scr_retval, __FILE__, __LINE__);
        }

        if (have_restart)
        {
            if (0 == rank)
            {
                printf("Restarting from %s...\n", dset);
            }

            #ifdef CKPT_DEBUG
                t1 = MPI_Wtime();
            #endif

            scr_retval = SCR_Start_restart(dset);

            #ifdef CKPT_DEBUG
                t2 = MPI_Wtime();
                t_start_restart += t2 - t1;
            #endif

            if (SCR_SUCCESS != scr_retval)
            {
                printf("%d: failed calling SCR_Start_restart: %d: @%s:%d\n",
                       rank, scr_retval, __FILE__, __LINE__);
            }

            //snprintf(path, sizeof(path), "%s/%s", dset, name);
            snprintf(path, sizeof(path), "%s/%s/%s", scr_prefix, dset, name);

            #ifdef CKPT_DEBUG
                t1 = MPI_Wtime();
            #endif

            scr_retval = SCR_Route_file(path, file);
            
            #ifdef CKPT_DEBUG
                t2 = MPI_Wtime();
                t_route_file_ch += t2 - t1;
            #endif

            if (SCR_SUCCESS != scr_retval)
            {
                printf("%d: failed calling SCR_Route_file: %d: @%s:%d\n",
                       rank, scr_retval, __FILE__, __LINE__);
            }
            else
            {
                if (read_ch(file, buf, length))
                {
                    found_checkpoint = 1;
                }
                else
                {
                    printf("%d: Could not read checkpoint %d from %s\n", rank, iteration, file);
                    found_checkpoint = 0;
                }
            }

            #ifdef CKPT_DEBUG
                t1 = MPI_Wtime();
            #endif

            scr_retval = SCR_Complete_restart(found_checkpoint);

            #ifdef CKPT_DEBUG
                t2 = MPI_Wtime();
                t_complete_restart += t2 - t1;
            #endif

            if (SCR_SUCCESS != scr_retval)
            {
                printf("%d: failed calling SCR_Complete_restart: %d: @%s:%d\n",
                       rank, scr_retval, __FILE__, __LINE__);
            }
            else
            {
                restarted = 1;
                ckpt_iteration = extract_final_number(dset);
                iteration = ckpt_iteration + 1;
                if (0 == rank)
                {
                    printf("Restarted from checkpoint %d\n", ckpt_iteration);
                }
            }
        }
    } while (have_restart && !restarted);
}

/**
 * Writes a checkpoint file.
 *
 * @param file   Pointer to the full path to a checkpoint file.
 * @param buf    Pointer to the input buffer.
 * @param length Quantity of elements in the input buffer.
 * @return       1 if the checkpoint file was written successfully; 0 otherwise.
 */
static int write_ch(char *file, TYPE *buf, int length)
{
    int rc, valid = 1;
    size_t return_value;
    FILE *pFile;

    /* open the file and write the checkpoint */
    pFile = fopen(file, "wb");
    if (NULL == pFile)
    {
        printf("%d: Could not open file %s\n", rank, file);
        valid = 0;
    }
    else
    {
        return_value = fwrite(buf, sizeof(TYPE), length, pFile);
        if (length != return_value)
        {
            valid = 0;
            printf("%d: Error writing %s\n", rank, file);
        }
        /* make sure the close is without error */
        rc = fclose(pFile);
        if (0 != rc)
        {
            valid = 0;
            printf("%d: Error closing %s\n", rank, file);
        }
    }
    return valid;
}

/**
 * Checks whether a checkpoint needs to be written and writes it.
 *
 * @param name   Pointer to the checkpoint filename.
 * @param buf    Pointer to the buffer.
 * @param length Quantity of elements in the input buffer.
 */
static void write_checkpoint(char *name, TYPE *buf, int length)
{
    int need_checkpoint = 0;
    int scr_retval, valid;
    char dirname[SCR_MAX_FILENAME];
    char path[SCR_MAX_FILENAME];
    char file[SCR_MAX_FILENAME];

    #ifdef CKPT_DEBUG
        t1 = MPI_Wtime();
    #endif

    #ifdef USE_SCR_NEED_CHECKPOINT
        scr_retval = SCR_Need_checkpoint(&need_checkpoint);

        #ifdef CKPT_DEBUG
            t2 = MPI_Wtime();
            t_need_checkpoint += t2 - t1;
        #endif

        if (SCR_SUCCESS != scr_retval)
        {
            printf("%d: failed calling SCR_Need_checkpoint: %d: @%s:%d\n",
                rank, scr_retval, __FILE__, __LINE__);
        }
    #else

        #ifdef CKPT_DEBUG
            t1 = MPI_Wtime();
        #endif

        need_checkpoint = (iteration % step == 0);

        #ifdef CKPT_DEBUG
            t2 = MPI_Wtime();
            t_need_checkpoint += t2 - t1;
        #endif

    #endif

    if (iteration == MAX_ITER-1) // last iteration
    {
        if (0 == rank) 
        {
            printf("Last iteration: will not save checkpoint\n");
        }
        need_checkpoint = 0;
    }

    if (need_checkpoint)
    {
        if (0 == rank)
        {
            printf("Writing checkpoint %d\n", iteration);
        }
        snprintf(dirname, sizeof(dirname), "timestep.%d", iteration);

        #ifdef CKPT_DEBUG
            t1 = MPI_Wtime();
        #endif

        scr_retval = SCR_Start_output(dirname, SCR_FLAG_CHECKPOINT);

        #ifdef CKPT_DEBUG
            t2 = MPI_Wtime();
            t_start_output += t2 - t1;
        #endif

        if (SCR_SUCCESS != scr_retval)
        {
            printf("%d: failed calling SCR_Start_output(): %d: @%s:%d\n",
                   rank, scr_retval, __FILE__, __LINE__);
        }

        snprintf(path, sizeof(path), "%s/%s/%s", scr_prefix, dirname, name);

        #ifdef CKPT_DEBUG
            t1 = MPI_Wtime();
        #endif

        #ifdef ONE_CKPT
            if (0 == rank)
            {
                scr_retval = SCR_Route_file(path, file);
            }
            else
            {
                scr_retval = SCR_SUCCESS;
            }
        #else
            scr_retval = SCR_Route_file(path, file);
        #endif

        #ifdef CKPT_DEBUG
            t2 = MPI_Wtime();
            t_route_file_out += t2 - t1;
        #endif

        if (SCR_SUCCESS != scr_retval)
        {
            printf("%d: failed calling SCR_Route_file(): %d: @%s:%d\n",
                rank, scr_retval, __FILE__, __LINE__);
        }

        #ifndef ONE_CKPT
            valid = write_ch(file, buf, length);
        #else
            if (0 == rank) 
            {
                valid = write_ch(file, buf, NB, MB);
            }
            else
            {
                valid = 1;
            }
        #endif

        #ifdef CKPT_DEBUG
            t1 = MPI_Wtime();
        #endif

        scr_retval = SCR_Complete_output(valid);

        #ifdef CKPT_DEBUG
            t2 = MPI_Wtime();
            t_complete_output += t2 - t1;
        #endif

        if (SCR_SUCCESS != scr_retval)
        {
            printf("%d: failed calling SCR_Complete_output: %d: @%s:%d\n",
                   rank, scr_retval, __FILE__, __LINE__);
        }
    }
}

/**
 * Prints the minimum and maximum timings of a specific loop in the program.
 *
 * @param scomm  MPI communicator for the processes involved in the timings.
 * @param rank   Rank of the current MPI process.
 * @param twf    Time (in seconds) taken for the specific loop in the current MPI process.
 */
void print_timings(MPI_Comm scomm, int rank, double twf)
{
    // Storage for min and max times
    double mtwf, Mtwf;

    // Perform reduction to find the minimum time across all MPI processes
    MPI_Reduce(&twf, &mtwf, 1, MPI_DOUBLE, MPI_MIN, 0, scomm);

    // Perform reduction to find the maximum time across all MPI processes
    MPI_Reduce(&twf, &Mtwf, 1, MPI_DOUBLE, MPI_MAX, 0, scomm);

    // If the current process is rank 0, print the min and max timings
    if (0 == rank)
    {
        printf(
            "##### Timings #####\n"
            "# MIN: %13.5e \t MAX: %13.5e\n",
            mtwf, Mtwf);
    }
}

/**
 * Performs one iteration of the Successive Over-Relaxation (SOR) method
 * on the input matrix and computes the squared L2-norm of the difference
 * between the input and output matrices.
 *
 * @param nm   Pointer to the output matrix after one iteration of the SOR method.
 * @param om   Pointer to the input matrix.
 * @param nb   Number of columns in the input matrix.
 * @param mb   Number of rows in the input matrix.
 * @return     The squared L2-norm of the difference between the input and output matrices.
 */
TYPE SOR1(TYPE *nm, TYPE *om, int nb, int mb)
{
    TYPE norm = 0.0;
    TYPE _W = 2.0 / (1.0 + M_PI / (TYPE)nb);
    int i, j, pos;

    // Iterate through each element of the matrix
    for (j = 0; j < mb; j++)
    {
        for (i = 0; i < nb; i++)
        {
            // Compute the position of the current element
            pos = 1 + i + (j + 1) * (nb + 2);

            // Update the current element using the SOR method
            nm[pos] = (1 - _W) * om[pos] +
                      _W / 4.0 * (nm[pos - 1] + om[pos + 1] + nm[pos - (nb + 2)] + om[pos + (nb + 2)]);

            // Accumulate the squared L2-norm of the difference
            norm += (nm[pos] - om[pos]) * (nm[pos] - om[pos]);
        }
    }

    return norm;
}

/**
 * Performs any required pre-initialization steps for the Jacobi method.
 * This function is a placeholder that can be extended if needed.
 *
 * @return     0 on successful completion.
 */
int preinit_jacobi_cpu(void)
{
    // Currently, there are no pre-initialization steps required for the
    // Jacobi method on the CPU. This function serves as a placeholder and
    // can be extended if necessary.

    return 0;
}

/**
 * Implements the Jacobi method for solving a system of linear equations using
 * MPI on a CPU. The convergence of the solution is controlled by the specified
 * epsilon value.
 *
 * @param matrix   Pointer to the input matrix of the linear system.
 * @param NB       Number of columns in the input matrix.
 * @param MB       Number of rows in the input matrix.
 * @param P        Number of partitions along the x-axis.
 * @param Q        Number of partitions along the y-axis.
 * @param comm     MPI communicator for the parallel computation.
 * @param epsilon  Convergence threshold for the Jacobi method.
 * @return         Number of iterations performed by the Jacobi method.
 */
int jacobi_cpu(TYPE *matrix, int NB, int MB, int P, int Q, MPI_Comm comm, TYPE epsilon)
{
    int ckpt_iteration, i;
    int size, ew_rank, ew_size, ns_rank, ns_size;
    TYPE *old_matrix, *new_matrix, *temp_matrix, *send_east, *send_west, *recv_east, *recv_west, diff_norm;
    double start_time, total_wf_time = 0; // timings
    char name[SCR_MAX_FILENAME];
    MPI_Comm ew, ns;

    MPI_Request req[8] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL,
                          MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};

    printf("Starting/resuming Jacobi method...\n");

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Query SCR informations from .scrconf (SCR_Config returns a pointer)
    // SCR_Finalize() will free the memory allocated for the SCR_Config variables
    {
        char *pstep;

        scr_prefix = (char *)SCR_Config("SCR_PREFIX");
        pstep = (char *)SCR_Config("SCR_CHECKPOINT_INTERVAL");
        step = atoi(pstep);
        free(pstep);
    }

    // Initialize SCR
    if (SCR_Init() != SCR_SUCCESS)
    {
        printf("SCR_Init failed\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
        return -1;
    }

    printf("Rank %d is joining the execution at iteration %d\n", rank, iteration);

    #ifndef ONE_CKPT
        snprintf(name, sizeof(name), "rank_%d.ckpt", rank);
    #else
        snprintf(name, sizeof(name), "checkpoint.ckpt");
    #endif

    old_matrix = matrix;
    new_matrix = (TYPE *)calloc(sizeof(TYPE), (NB + 2) * (MB + 2));
    send_east = (TYPE *)malloc(sizeof(TYPE) * MB);
    send_west = (TYPE *)malloc(sizeof(TYPE) * MB);
    recv_east = (TYPE *)malloc(sizeof(TYPE) * MB);
    recv_west = (TYPE *)malloc(sizeof(TYPE) * MB);

    // Create the north-south and east-west communicators
    MPI_Comm_split(comm, rank % P, rank, &ns);
    MPI_Comm_size(ns, &ns_size);
    MPI_Comm_rank(ns, &ns_rank);

    MPI_Comm_split(comm, rank / P, rank, &ew);
    MPI_Comm_size(ew, &ew_size);
    MPI_Comm_rank(ew, &ew_rank);

    try_restart(name, &old_matrix, (NB + 2) * (MB + 2));

    start_time = MPI_Wtime();
    do
    {
        // Post receives from the neighbors
        if (0 != ns_rank)
            MPI_Irecv(RECV_NORTH(old_matrix), NB, MPI_TYPE, ns_rank - 1, 0, ns, &req[0]);
        if ((ns_size - 1) != ns_rank)
            MPI_Irecv(RECV_SOUTH(old_matrix), NB, MPI_TYPE, ns_rank + 1, 0, ns, &req[1]);
        if ((ew_size - 1) != ew_rank)
            MPI_Irecv(recv_east, MB, MPI_TYPE, ew_rank + 1, 0, ew, &req[2]);
        if (0 != ew_rank)
            MPI_Irecv(recv_west, MB, MPI_TYPE, ew_rank - 1, 0, ew, &req[3]);

        // Post the sends
        if (0 != ns_rank)
            MPI_Isend(SEND_NORTH(old_matrix), NB, MPI_TYPE, ns_rank - 1, 0, ns, &req[4]);
        if ((ns_size - 1) != ns_rank)
            MPI_Isend(SEND_SOUTH(old_matrix), NB, MPI_TYPE, ns_rank + 1, 0, ns, &req[5]);
        for (i = 0; i < MB; i++)
        {
            send_west[i] = old_matrix[(i + 1) * (NB + 2) + 1];      // The real local data
            send_east[i] = old_matrix[(i + 1) * (NB + 2) + NB + 0]; // Not the ghost region
        }
        if ((ew_size - 1) != ew_rank)
            MPI_Isend(send_east, MB, MPI_TYPE, ew_rank + 1, 0, ew, &req[6]);
        if (0 != ew_rank)
            MPI_Isend(send_west, MB, MPI_TYPE, ew_rank - 1, 0, ew, &req[7]);

        // Wait until they all complete
        MPI_Waitall(8, req, MPI_STATUSES_IGNORE);
        for (i = 0; i < MB; i++)
        {
            old_matrix[(i + 1) * (NB + 2)] = recv_west[i];
            old_matrix[(i + 1) * (NB + 2) + NB + 1] = recv_east[i];
        }

        // Replicate the east-west newly received data
        for (i = 0; i < MB; i++)
        {
            new_matrix[(i + 1) * (NB + 2)] = old_matrix[(i + 1) * (NB + 2)];
            new_matrix[(i + 1) * (NB + 2) + NB + 1] = old_matrix[(i + 1) * (NB + 2) + NB + 1];
        }

        // Replicate the north-south neighbors
        for (i = 0; i < NB; i++)
        {
            new_matrix[i + 1] = old_matrix[i + 1];
            new_matrix[(NB + 2) * (MB + 1) + i + 1] = old_matrix[(NB + 2) * (MB + 1) + i + 1];
        }

        diff_norm = SOR1(new_matrix, old_matrix, NB, MB);

        if (verbose)
            printf("Rank %d norm %f at iteration %d\n", rank, diff_norm, iteration);

        // Allreduce to compute the global norm
        MPI_Allreduce(MPI_IN_PLACE, &diff_norm, 1, MPI_TYPE, MPI_SUM, comm);

        if (0 == rank)
        {
            printf("Iteration %4d norm %f\n", iteration, sqrtf(diff_norm));
        }

        // Swap the old and new matrices
        temp_matrix = old_matrix;
        old_matrix = new_matrix;
        new_matrix = temp_matrix;

        write_checkpoint(name, old_matrix, (NB + 2) * (MB + 2));

        // Increment the iteration
        iteration++;

    } while ((iteration < MAX_ITER) && (sqrt(diff_norm) > epsilon));

    total_wf_time = MPI_Wtime() - start_time;

    print_timings(comm, rank, total_wf_time);

    #ifdef CKPT_DEBUG
        double avg_t_have_restart;
        MPI_Reduce(&t_have_restart, &avg_t_have_restart, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        avg_t_have_restart /= size;

        double avg_t_start_restart;
        MPI_Reduce(&t_start_restart, &avg_t_start_restart, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        avg_t_start_restart /= size;

        double avg_t_route_file_ch;
        MPI_Reduce(&t_route_file_ch, &avg_t_route_file_ch, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        avg_t_route_file_ch /= size;

        double avg_t_complete_restart;
        MPI_Reduce(&t_complete_restart, &avg_t_complete_restart, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        avg_t_complete_restart /= size;

        double avg_t_need_checkpoint;
        MPI_Reduce(&t_need_checkpoint, &avg_t_need_checkpoint, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        avg_t_need_checkpoint /= size;

        double avg_t_start_output;
        MPI_Reduce(&t_start_output, &avg_t_start_output, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        avg_t_start_output /= size;

        double avg_t_route_file_out;
        MPI_Reduce(&t_route_file_out, &avg_t_route_file_out, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        avg_t_route_file_out /= size;

        double avg_t_complete_output;
        MPI_Reduce(&t_complete_output, &avg_t_complete_output, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        avg_t_complete_output /= size;
        
        // If the current process is rank 0, print the avg timings
        if (0 == rank)
        {
            printf("##### Debug timings #####\n");
            printf("# t_have_restart     (AVG): %13.5e\n", avg_t_have_restart);
            printf("# t_start_restart    (AVG): %13.5e\n", avg_t_start_restart);
            printf("# t_route_file_ch    (AVG): %13.5e\n", avg_t_route_file_ch);
            printf("# t_complete_restart (AVG): %13.5e\n", avg_t_complete_restart);
        #ifdef USE_SCR_NEED_CHECKPOINT
            printf("# t_need_checkpoint  (AVG): %13.5e\n", avg_t_need_checkpoint);
        #else
            printf("# manual_ch_check    (AVG): %13.5e\n", avg_t_need_checkpoint);
        #endif
            printf("# t_start_output     (AVG): %13.5e\n", avg_t_start_output);
            printf("# t_route_file_out   (AVG): %13.5e\n", avg_t_route_file_out);
            printf("# t_complete_output  (AVG): %13.5e\n", avg_t_complete_output);
        }
    #endif

    // Free the memory allocated for matrices and buffers
    // If the 'matrix' variable is different from 'old_matrix', free 'old_matrix'; otherwise, free 'new_matrix'
    free(matrix != old_matrix ? old_matrix : new_matrix);
    
    // Free the memory allocated for send and receive buffers
    free(send_west);
    free(send_east);
    free(recv_west);
    free(recv_east);

    MPI_Comm_free(&ns);
    MPI_Comm_free(&ew);

    SCR_Finalize();

    return iteration;
}
