#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>

#include "../inc/argument_utils.h"

#include <mpi.h>

#define WALLTIME(t) ((double)(t).tv_sec + 1e-6 * (double)(t).tv_usec)

MPI_Comm
    cart; // holds the cartesian communicator
MPI_Datatype
    grid,
    subgrid,
    y_column,
    x_column;

int
    rank,
    comm_size,
    local_rows,
    local_cols;

// Cartesian information
int n_dims = 2; // number of dimensions for cartesian communicator
int periods[2] = {0, 0}; // periodic information for dimensions
int *dims; // pointer to array of number of processes in each dimension

#define MPI_RANK_ROOT  ( rank == 0 )

struct timeval
    t_start,
    t_stop;
double
    t_total;

typedef int64_t int_t;
typedef double real_t;

int_t
    N,
    max_iteration,
    snapshot_frequency;

const real_t
    domain_size= 10.0,
    gravity = 9.81,
    density = 997.0;

real_t
    *mass[2] = { NULL, NULL},
    *mass_velocity_x[2] = { NULL, NULL },
    *mass_velocity_y[2] = { NULL, NULL },
    *mass_velocity = NULL,
    *velocity_x = NULL,
    *velocity_y = NULL,
    *acceleration_x = NULL,
    *acceleration_y = NULL,
    dx,
    dt;

#define PN(y,x)         mass[0][(y)*(local_cols+2)+(x)] // Mass at time t
#define PN_next(y,x)    mass[1][(y)*(local_cols+2)+(x)] // Mass at time t+1
#define PNU(y,x)        mass_velocity_x[0][(y)*(local_cols+2)+(x)] // Mass velocity in x-drection at time t
#define PNU_next(y,x)   mass_velocity_x[1][(y)*(local_cols+2)+(x)] // Mass velocity in x-drection at time t + 1
#define PNV(y,x)        mass_velocity_y[0][(y)*(local_cols+2)+(x)]
#define PNV_next(y,x)   mass_velocity_y[1][(y)*(local_cols+2)+(x)]
#define PNUV(y,x)       mass_velocity[(y)*(local_cols+2)+(x)] 
#define U(y,x)          velocity_x[(y)*(local_cols+2)+(x)]
#define V(y,x)          velocity_y[(y)*(local_cols+2)+(x)]
#define DU(y,x)         acceleration_x[(y)*(local_cols+2)+(x)]
#define DV(y,x)         acceleration_y[(y)*(local_cols+2)+(x)]

void time_step ( void );
void boundary_condition ( real_t *domain_variable, int sign );
void create_types ( void );
void domain_init ( void );
void domain_save ( int_t iteration );
void domain_finalize ( void );


void
swap ( real_t** t1, real_t** t2 )
{
    real_t* tmp;
	tmp = *t1;
	*t1 = *t2;
	*t2 = tmp;
}


int
main ( int argc, char **argv )
{
    MPI_Init ( &argc, &argv );

    // TODO 1 Create a communicator with cartesian topology
    MPI_Comm_size ( MPI_COMM_WORLD, &comm_size );
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );  
    
    // Allocate processes to grid
    dims = calloc(n_dims, sizeof(int)); // allocate adequately sized array   
    MPI_Dims_create(comm_size, n_dims, dims); // find number of processes in each dimension
    MPI_Cart_create(MPI_COMM_WORLD, n_dims, dims, periods, 0, &cart); 

    if ( MPI_RANK_ROOT )
    {
        OPTIONS *options = parse_args( argc, argv );
        if ( !options )
        {
            fprintf( stderr, "Argument parsing failed\n" );
            exit(1);
        }

        N = options->N;
        max_iteration = options->max_iteration;
        snapshot_frequency = options->snapshot_frequency;

    }
    

    MPI_Bcast ( &N, 1, MPI_INT64_T, 0, MPI_COMM_WORLD );
    MPI_Bcast ( &max_iteration, 1, MPI_INT64_T, 0, MPI_COMM_WORLD );
    MPI_Bcast ( &snapshot_frequency, 1, MPI_INT64_T, 0, MPI_COMM_WORLD );
    

    // TODO 2 Find the number of columns and rows of each subgrid
    //        and find the local x and y offsets for each process' subgrid
    // Uncomment domain_save() and create_types() after TODO2 is complete
    domain_init();

    create_types();

    gettimeofday ( &t_start, NULL );

    for ( int_t iteration = 0; iteration<=max_iteration; iteration++ )
    {
        // TODO 5 Implement border exchange
        int 
            y_neighbor_up, y_neighbor_down,
            x_neighbor_left, x_neighbor_right;
        
        MPI_Cart_shift(cart, 0, 1, &y_neighbor_up, &y_neighbor_down);
        MPI_Cart_shift(cart, 1, 1, &x_neighbor_left, &x_neighbor_right);

        if (y_neighbor_down != MPI_PROC_NULL) { 
            MPI_Sendrecv(&DU(1, 1), 1, x_column, y_neighbor_down, 0, 
                         &DU(0, 1), 1, x_column, y_neighbor_down, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&DV(1, 1), 1, x_column, y_neighbor_down, 0, 
                         &DU(0, 1), 1, x_column, y_neighbor_down, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&PN(1, 1), 1, x_column, y_neighbor_down, 0, 
                         &PN(0, 1), 1, x_column, y_neighbor_down, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&PNU(1, 1), 1, x_column, y_neighbor_down, 0, 
                         &PNU(0, 1), 1, x_column, y_neighbor_down, 0, 
                         cart, MPI_STATUS_IGNORE);             
            MPI_Sendrecv(&PNV(1, 1), 1, x_column, y_neighbor_down, 0, 
                         &PNV(0, 1), 1, x_column, y_neighbor_down, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&PNUV(1, 1), 1, x_column, y_neighbor_down, 0, 
                         &PNUV(0, 1), 1, x_column, y_neighbor_down, 0, 
                         cart, MPI_STATUS_IGNORE);
        } else {
            MPI_Sendrecv(&DU(local_rows - 1, 1), 1, x_column, y_neighbor_up, 0, 
                         &DU(local_rows + 1, 1), 1, x_column, y_neighbor_up, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&DV(local_rows - 1, 1), 1, x_column, y_neighbor_up, 0, 
                         &DV(local_rows + 1, 1), 1, x_column, y_neighbor_up, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&PN(local_rows - 1, 1), 1, x_column, y_neighbor_up, 0, 
                         &PN(local_rows + 1, 1), 1, x_column, y_neighbor_up, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&PNU(local_rows - 1, 1), 1, x_column, y_neighbor_up, 0, 
                         &PNU(local_rows + 1, 1), 1, x_column, y_neighbor_up, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&PNV(local_rows - 1, 1), 1, x_column, y_neighbor_up, 0, 
                         &PNV(local_rows + 1, 1), 1, x_column, y_neighbor_up, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&PNUV(local_rows - 1, 1), 1, x_column, y_neighbor_up, 0, 
                         &PNUV(local_rows + 1, 1), 1, x_column, y_neighbor_up, 0, 
                         cart, MPI_STATUS_IGNORE);
        }

        if (x_neighbor_right != MPI_PROC_NULL) {
            MPI_Sendrecv(&DU(1, local_cols), 1, y_column, x_neighbor_right, 0, 
                         &DU(1, local_cols + 1), 1, y_column, x_neighbor_right, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&DV(1, local_cols), 1, y_column, x_neighbor_right, 0, 
                         &DV(1, local_cols + 1), 1, y_column, x_neighbor_right, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&PN(1, local_cols), 1, y_column, x_neighbor_right, 0, 
                         &PN(1, local_cols + 1), 1, y_column, x_neighbor_right, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&PNU(1, local_cols), 1, y_column, x_neighbor_right, 0, 
                         &PNU(1, local_cols + 1), 1, y_column, x_neighbor_right, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&PNV(1, local_cols), 1, y_column, x_neighbor_right, 0, 
                         &PNV(1, local_cols + 1), 1, y_column, x_neighbor_right, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&PNUV(1, local_cols), 1, y_column, x_neighbor_right, 0, 
                         &PNUV(1, local_cols + 1), 1, y_column, x_neighbor_right, 0, 
                         cart, MPI_STATUS_IGNORE);
        } else {
            MPI_Sendrecv(&DU(1, 1), 1, y_column, x_neighbor_left, 0, 
                         &DU(1, 0), 1, y_column, x_neighbor_left, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&DV(1, 1), 1, y_column, x_neighbor_left, 0, 
                         &DV(1, 0), 1, y_column, x_neighbor_left, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&PN(1, 1), 1, y_column, x_neighbor_left, 0, 
                         &PN(1, 0), 1, y_column, x_neighbor_left, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&PNU(1, 1), 1, y_column, x_neighbor_left, 0, 
                         &PNU(1, 0), 1, y_column, x_neighbor_left, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&PNV(1, 1), 1, y_column, x_neighbor_left, 0, 
                         &PNV(1, 0), 1, y_column, x_neighbor_left, 0, 
                         cart, MPI_STATUS_IGNORE);
            MPI_Sendrecv(&PNUV(1, 1), 1, y_column, x_neighbor_left, 0, 
                         &PNUV(1, 0), 1, y_column, x_neighbor_left, 0, 
                         cart, MPI_STATUS_IGNORE);
        }

        // TODO 4 Change application of boundary condition to match cartesian topology
        boundary_condition ( mass[0], 1 );
        boundary_condition ( mass_velocity_x[0], -1 );
        boundary_condition ( mass_velocity_y[0], -1 );

        // TODO 3 Update the area of iteration in the time step
        time_step();

        if ( iteration % snapshot_frequency == 0 )
        {
            if ( MPI_RANK_ROOT )
            {
                printf (
                    "Iteration %lld of %lld, (%.2lf%% complete)\n",
                    iteration,
                    max_iteration,
                    100.0 * (real_t) iteration / (real_t) max_iteration
                );

            }   
            domain_save ( iteration );
        }

        swap ( &mass[0], &mass[1] );
        swap ( &mass_velocity_x[0], &mass_velocity_x[1] );
        swap ( &mass_velocity_y[0], &mass_velocity_y[1] );
    }

    domain_finalize();

    gettimeofday ( &t_stop, NULL );
    t_total = WALLTIME(t_stop) - WALLTIME(t_start);

    if ( MPI_RANK_ROOT )
        printf ( "%.2lf seconds total runtime\n", t_total );

    MPI_Finalize();

    exit ( EXIT_SUCCESS );
}


void
time_step ( void )
{
    // TODO 3 Update the area of iteration in the time step

    for ( int_t y=1; y<=local_rows; y++ )
        for ( int_t x=1; x<=local_cols; x++ )
        {
            U(y,x) = PNU(y,x) / PN(y,x);
            V(y,x) = PNV(y,x) / PN(y,x);
        }

    for ( int_t y=1; y<=local_rows; y++ )
        for ( int_t x=1; x<=local_cols; x++ )
        {
            PNUV(y,x) = PN(y,x) * U(y,x) * V(y,x);
        }

    for ( int_t y=0; y<=local_rows+1; y++ )
        for ( int_t x=0; x<=local_cols+1; x++ )
        {
            DU(y,x) = PN(y,x) * U(y,x) * U(y,x)
                    + 0.5 * gravity * ( PN(y,x) * PN(y,x) / density );
            DV(y,x) = PN(y,x) * V(y,x) * V(y,x)
                    + 0.5 * gravity * ( PN(y,x) * PN(y,x) / density );
        }

    for ( int_t y=1; y<=local_rows; y++ )
        for ( int_t x=1; x<=local_cols; x++ )
        {
            PNU_next(y,x) = 0.5*( PNU(y,x+1) + PNU(y,x-1) ) - dt*(
                            ( DU(y,x+1) - DU(y,x-1) ) / (2*dx)
                          + ( PNUV(y,x+1) - PNUV(y,x-1) ) / (2*dx)
            );
        }

    for ( int_t y=1; y<=local_rows; y++ )
        for ( int_t x=1; x<=local_cols; x++ )
        {
            PNV_next(y,x) = 0.5*( PNV(y+1,x) + PNV(y-1,x) ) - dt*(
                            ( DV(y+1,x) - DV(y-1,x) ) / (2*dx)
                          + ( PNUV(y+1,x) - PNUV(y-1,x) ) / (2*dx)
            );
        }

    for ( int_t y=1; y<=local_rows; y++ )
        for ( int_t x=1; x<=local_cols; x++ )
        {
            PN_next(y,x) = 0.25*( PN(y,x+1) + PN(y,x-1) + PN(y+1,x) + PN(y-1,x) ) - dt*(
                           ( PNU(y,x+1) - PNU(y,x-1) ) / (2*dx)
                         + ( PNV(y+1,x) - PNV(y-1,x) ) / (2*dx)
            );
        }
}


void
boundary_condition ( real_t *domain_variable, int sign )
{
    // TODO 4 Change application of boundary condition to match cartesian topology

    #define VAR(y,x) domain_variable[(y)*(local_cols+2)+(x)]

    VAR(   0, 0   ) = sign*VAR(   2, 2   );
    VAR( local_rows+1, 0   ) = sign*VAR( local_rows-1, 2   );
    VAR(   0, local_cols+1 ) = sign*VAR(   2, local_cols-1 );
    VAR( local_rows+1, local_cols+1 ) = sign*VAR( local_rows-1, local_cols-1 );

    for ( int_t y=1; y<=local_rows; y++ ) VAR(   y, 0   ) = sign*VAR(   y, 2   );
    for ( int_t y=1; y<=local_rows; y++ ) VAR(   y, local_cols+1 ) = sign*VAR(   y, local_cols-1 );
    for ( int_t x=1; x<=local_cols; x++ ) VAR(   0, x   ) = sign*VAR(   2, x   );
    for ( int_t x=1; x<=local_cols; x++ ) VAR( local_rows+1, x   ) = sign*VAR( local_rows-1, x   );
    
    #undef VAR
}


void
create_types ( void )
{
    int cart_rank, cart_offset[2];
    MPI_Comm_rank ( cart, &cart_rank );
    MPI_Cart_coords ( cart, cart_rank, 2, cart_offset);

    MPI_Type_create_subarray ( 2, (int[2]) { local_rows+2, local_cols+2 }, (int[2]) { local_rows, local_cols }, (int[2]) {1,1}, MPI_ORDER_C, MPI_DOUBLE, &subgrid );
    MPI_Type_create_subarray ( 2, (int[2]) {N, N} , (int[2]) { local_rows, local_cols }, (int[2]) { cart_offset[0] * local_rows, cart_offset[1] * local_cols}, MPI_ORDER_C, MPI_DOUBLE, &grid );

    MPI_Type_commit ( &subgrid );
    MPI_Type_commit ( &grid ) ;

    // create column type
    MPI_Type_vector(local_rows, 1, local_cols + 2, MPI_DOUBLE, &y_column);
    MPI_Type_vector(local_cols, 1, 1, MPI_DOUBLE, &x_column);

    MPI_Type_commit(&y_column);
    MPI_Type_commit(&x_column);

}


void
domain_init ( void )
{
    // TODO 2 Find the number of columns and rows of each subgrid
    // Hint: you can get useful information from the cartesian communicator
    local_rows  = N / dims[0]; // number of rows = discretization steps / number of processes in dim 0 (evenly divisible) 
    local_cols  = N / dims[1]; // number of cols = discretization steps / number of processes in dim 1 (evenly divisible)

    int_t local_size = (local_rows + 2) * (local_cols + 2);

    mass[0] = calloc ( local_size, sizeof(real_t) );
    mass[1] = calloc ( local_size, sizeof(real_t) );

    mass_velocity_x[0] = calloc ( local_size, sizeof(real_t) );
    mass_velocity_x[1] = calloc ( local_size, sizeof(real_t) );
    mass_velocity_y[0] = calloc ( local_size, sizeof(real_t) );
    mass_velocity_y[1] = calloc ( local_size, sizeof(real_t) );

    mass_velocity = calloc ( local_size, sizeof(real_t) );

    velocity_x = calloc ( local_size, sizeof(real_t) );
    velocity_y = calloc ( local_size, sizeof(real_t) );

    acceleration_x = calloc ( local_size, sizeof(real_t) );
    acceleration_y = calloc ( local_size, sizeof(real_t) );

    // TODO 2 Find the local x and y offsets for each process' subgrid
    // Hint: you can get useful information from the cartesian communicator
    int coords[n_dims]; // store process coordinates

    MPI_Cart_coords(cart, rank, n_dims, coords);

    int_t local_x_offset = coords[1] * (N / dims[1]); // x-coord * N
    int_t local_y_offset = coords[0] * ((N / dims[0])); // y-coord * N
    

    for ( int_t y=1; y<=local_rows; y++ )
    {
        for ( int_t x=1; x<=local_cols; x++ )
        {
            PN(y,x) = 1e-3;
            PNU(y,x) = 0.0;
            PNV(y,x) = 0.0;

            real_t cx = (local_x_offset + x) - N/2;
            real_t cy = (local_y_offset + y) - N/2;

            if ( sqrt ( cx*cx + cy*cy ) < N/20.0 )
            {
                PN(y,x) -= 5e-4*exp (
                    - 4*pow( cx, 2.0 ) / (real_t)(N)
                    - 4*pow( cy, 2.0 ) / (real_t)(N)
                );
            }

            PN(y,x) *= density;
        }
    }

    dx = domain_size / (real_t) N;
    dt = 5e-2;

    free(dims); // no longer needed
}


void
domain_save ( int_t iteration )
{
    int_t index = iteration / snapshot_frequency;
    char filename[256];
    memset ( filename, 0, 256*sizeof(char) );
    sprintf ( filename, "data/%.5lld.bin", index );

    MPI_File out;
    MPI_File_open ( cart, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &out );

    MPI_File_set_view ( out, 0, MPI_DOUBLE, grid, "native", MPI_INFO_NULL );
    MPI_File_write_all ( out, mass[0], 1, subgrid, MPI_STATUS_IGNORE );

    MPI_File_close ( &out );
}


void
domain_finalize ( void )
{
    free ( mass[0] );
    free ( mass[1] );
    free ( mass_velocity_x[0] );
    free ( mass_velocity_x[1] );
    free ( mass_velocity_y[0] );
    free ( mass_velocity_y[1] );
    free ( mass_velocity );
    free ( velocity_x );
    free ( velocity_y );
    free ( acceleration_x );
    free ( acceleration_y );
}
