#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>

#include "../inc/argument_utils.h"

typedef int64_t int_t;
typedef double real_t;

int_t
    N,
    max_iteration,
    snapshot_frequency;

const real_t
    domain_size = 10.0,
    gravity = 9.81,
    density = 997.0;

// Global MPI variables
int size, rank, grid_size;

// Grid variables
int offset, 
    east_neighbor, 
    west_neighbor;


real_t
    *mass[2] = { NULL, NULL },
    *mass_velocity_x[2] = { NULL, NULL },
    *velocity_x = NULL,
    *acceleration_x = NULL,
    dx,
    dt;

#define PN(x)        mass[0][(x)]
#define PN_next(x)   mass[1][(x)]
#define PNU(x)       mass_velocity_x[0][(x)]
#define PNU_next(x)  mass_velocity_x[1][(x)]
#define U(x)         velocity_x[(x)]
#define DU(x)        acceleration_x[(x)]

void time_step ( void );
void boundary_condition( real_t *domain_variable, int sign, int rank );
void domain_init ( int rank );
void domain_save ( int_t iteration, int rank );
void domain_finalize ( void );


void
swap ( real_t** m1, real_t** m2 )
{
    real_t* tmp;
    tmp = *m1;
    *m1 = *m2;
    *m2 = tmp;
}


int
main ( int argc, char **argv )
{
    // TODO 1 Initialize MPI

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // TODO 2 Parse arguments in the rank 0 processes
    // and broadcast to other processes

    OPTIONS *options;

    if (rank == 0) {
        options = parse_args( argc, argv );

        if ( !options )
        {
            fprintf( stderr, "Argument parsing failed\n" );
            exit(1);
        }

        N = options->N;
        max_iteration = options->max_iteration;
        snapshot_frequency = options->snapshot_frequency;

    }
    
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_iteration, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&snapshot_frequency, 1, MPI_INT, 0, MPI_COMM_WORLD);    
    
    grid_size = N / size;
    
    east_neighbor = rank + 1;
    west_neighbor = rank - 1;

    // we separate the if-checks in case of single process
    if (rank == 0) {
        west_neighbor = size - 1;
    } 
    
    if (rank == size - 1) {
        east_neighbor = 0;
    }

    // TODO 3 Allocate space for each process' sub-grid
    // and initialize data for the sub-grid

    domain_init(rank);

    for ( int_t iteration = 0; iteration <= max_iteration; iteration++ )
    {
        // TODO 7 Communicate border values 
        
        // east border values
        MPI_Send(&PN(grid_size), 1, MPI_DOUBLE, east_neighbor, 0, MPI_COMM_WORLD);
        MPI_Recv(&PN(0), 1, MPI_DOUBLE, west_neighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        MPI_Send(&PNU(grid_size), 1, MPI_DOUBLE, east_neighbor, 0, MPI_COMM_WORLD); 
        MPI_Recv(&PNU(0), 1, MPI_DOUBLE, west_neighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


        // send velocity if not on east border
        if (rank != size - 1) {
            MPI_Send(&U(grid_size), 1, MPI_DOUBLE, east_neighbor, 0, MPI_COMM_WORLD);
        }

        // receive velocity if not on west border
        if (rank != 0) {
            MPI_Recv(&U(0), 1, MPI_DOUBLE, west_neighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // west border values
        MPI_Send(&PN(1), 1, MPI_DOUBLE, west_neighbor, 0, MPI_COMM_WORLD);
        MPI_Recv(&PN(grid_size + 1), 1, MPI_DOUBLE, east_neighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        MPI_Send(&PNU(1), 1, MPI_DOUBLE, west_neighbor, 0, MPI_COMM_WORLD);
        MPI_Recv(&PNU(grid_size + 1), 1, MPI_DOUBLE, east_neighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (rank != 0) {
            MPI_Send(&U(1), 1, MPI_DOUBLE, west_neighbor, 0, MPI_COMM_WORLD);
        }

        if (rank != size - 1) {
            MPI_Recv(&U(grid_size + 1), 1, MPI_DOUBLE, east_neighbor, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    
        // TODO 5 Boundary conditions
        boundary_condition(mass[0], 1, rank);
        boundary_condition(mass_velocity_x[0], -1, rank);
        
        // TODO 4 Time step calculations
        time_step();



        if ( iteration % snapshot_frequency == 0 )
        {
            printf (
                "[RANK %d] Iteration %ld of %ld (%.2lf%% complete)\n",
                rank,
                iteration,
                max_iteration,
                100.0 * (real_t) iteration / (real_t) max_iteration
            );

            // TODO 6 MPI I/O
            domain_save ( iteration, rank);
        }

        swap( &mass[0], &mass[1] );
        swap( &mass_velocity_x[0], &mass_velocity_x[1] );
    }

    domain_finalize();

    // TODO 1 Finalize MPI

    MPI_Finalize();

    exit ( EXIT_SUCCESS );
}


void
time_step ( void )
{
    // TODO 4 Time step calculations


    for ( int_t x=0; x<=grid_size + 1; x++ )
    {
        DU(x) = PN(x) * U(x) * U(x) + 0.5 * gravity * PN(x) * PN(x) / density;
        /**
        if (x == 0 || x == grid_size + 1) {
            printf("[RANK %d] Current acceleration at position %d: %lld\n", rank, x, DU(x));
            printf("[RANK %d] Current velocity at position %d: %lld\n", rank, x, U(x));
            printf("[RANK %d] Current mass at position %d: %lld\n", rank, x, PN(x));

        }
        */
    }

    for ( int_t x=1; x<=grid_size; x++ )
    {
        PNU_next(x) = 0.5*( PNU(x+1) + PNU(x-1) ) - dt*(
                      ( DU(x+1) - DU(x-1) ) / (2*dx)
        );
    }

    for ( int_t x=1; x<=grid_size; x++ )
    {
        PN_next(x) = 0.5*( PN(x+1) + PN(x-1) ) - dt*(
                       ( PNU(x+1) - PNU(x-1) ) / (2*dx)
        );
    }

    for ( int_t x=1; x<=grid_size; x++ )
    {
        U(x) = PNU_next(x) / PN_next(x);
    }
}


void
boundary_condition ( real_t *domain_variable, int sign, int rank )
{
    // TODO 5 Boundary conditions

    #define VAR(x) domain_variable[(x)]
    if (rank == 0) {
        VAR( 0 ) = sign*VAR( 2 );
    } 
    
    if (rank == size - 1) {
        VAR( grid_size + 1 ) = sign*VAR( grid_size-1 );
    };
    #undef VAR
}


void
domain_init ( int rank )
{
    // TODO 3 Allocate space for each process' sub-grid
    // and initialize data for the sub-grid

    mass[0] = calloc ( (grid_size+2), sizeof(real_t) ); // PN(x)
    mass[1] = calloc ( (grid_size+2),  sizeof(real_t) ); // PN_next(x)

    mass_velocity_x[0] = calloc ( (grid_size+2), sizeof(real_t) ); // PNU(x)
    mass_velocity_x[1] = calloc ( (grid_size+2),  sizeof(real_t) ); // PNU_next(x)

    velocity_x = calloc ( (grid_size+2), sizeof(real_t) );
    acceleration_x = calloc ( (grid_size+2), sizeof(real_t) );

    // each process is allocated a grid of responsibility based on their rank
    // 0 is first, 3 is last
    offset = grid_size * rank;  

    // Data initialization
    for ( int_t x=1; x<=grid_size; x++ )
    {
        PN(x) = 1e-3;
        PNU(x) = 0.0;

        real_t c = x + offset - N/2; // add the offset to check the position in the global grid
        if ( sqrt ( c*c ) < N/20.0 )
        {
            PN(x) -= 5e-4*exp (
                    - 4*pow( c, 2.0 ) / (real_t)(N)
            );
        }

        PN(x) *= density;
    }

    dx = domain_size / (real_t) N;
    dt = 0.1*dx;

}


void
domain_save ( int_t iteration, int rank )
{
    int_t index = iteration / snapshot_frequency;
    char filename[256];
    memset ( filename, 0, 256*sizeof(char) );
    sprintf ( filename, "data/%.5ld.bin", index );

    // TODO 6 MPI I/O

    MPI_File out;

    MPI_File_open(MPI_COMM_WORLD,
                  filename,
                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL,
                  &out);

    if ( ! out ) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        exit(1);
    }
    

    MPI_Offset offset = rank * grid_size * sizeof(real_t);
    MPI_File_write_at_all(out, offset, &mass[0][1], grid_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&out);

}   


void
domain_finalize ( void )
{
    free ( mass[0] );
    free ( mass[1] );
    free ( mass_velocity_x[0] );
    free ( mass_velocity_x[1] );
    free ( velocity_x );
    free ( acceleration_x );
}

