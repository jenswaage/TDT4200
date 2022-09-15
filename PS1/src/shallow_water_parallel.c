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
void boundary_condition( real_t *domain_variable, int sign );
void domain_init ( void );
void domain_save ( int_t iteration );
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
    int size, rank;
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
        
    }

    MPI_Bcast(options, sizeof(OPTIONS), MPI_INT, 0, MPI_COMM_WORLD);

    N = options->N;
    max_iteration = options->max_iteration;
    snapshot_frequency = options->snapshot_frequency;

    printf("[RANK %d] {N: %d, max_iteration: %d, snapshot_frequency: %d}\n", rank, N, max_iteration, snapshot_frequency);
    // TODO 3 Allocate space for each process' sub-grid
    // and initialize data for the sub-grid
    domain_init();

    for ( int_t iteration = 0; iteration <= max_iteration; iteration++ )
    {
        // TODO 7 Communicate border values


        // TODO 5 Boundary conditions
        boundary_condition(mass[0], 1);
        boundary_condition(mass_velocity_x[0], -1);

        // TODO 4 Time step calculations
        time_step();

        if ( iteration % snapshot_frequency == 0 )
        {
            printf (
                "Iteration %ld of %ld (%.2lf%% complete)\n",
                iteration,
                max_iteration,
                100.0 * (real_t) iteration / (real_t) max_iteration
            );

            // TODO 6 MPI I/O
            domain_save ( iteration );
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

    for ( int_t x=0; x<=N+1; x++ )
    {
        DU(x) = PN(x) * U(x) * U(x)
                + 0.5 * gravity * PN(x) * PN(x) / density;
    }

    for ( int_t x=1; x<=N; x++ )
    {
        PNU_next(x) = 0.5*( PNU(x+1) + PNU(x-1) ) - dt*(
                      ( DU(x+1) - DU(x-1) ) / (2*dx)
        );
    }

    for ( int_t x=1; x<=N; x++ )
    {
        PN_next(x) = 0.5*( PN(x+1) + PN(x-1) ) - dt*(
                       ( PNU(x+1) - PNU(x-1) ) / (2*dx)
        );
    }

    for ( int_t x=1; x<=N; x++ )
    {
        U(x) = PNU_next(x) / PN_next(x);
    }
}


void
boundary_condition ( real_t *domain_variable, int sign )
{
    // TODO 5 Boundary conditions

    #define VAR(x) domain_variable[(x)]
    VAR(   0 ) = sign*VAR( 2   );
    VAR( N+1 ) = sign*VAR( N-1 );
    #undef VAR
}


void
domain_init ( void )
{
    // TODO 3 Allocate space for each process' sub-grid
    // and initialize data for the sub-grid

    mass[0] = calloc ( (N+2), sizeof(real_t) );
    mass[1] = calloc ( (N+2),  sizeof(real_t) );

    mass_velocity_x[0] = calloc ( (N+2), sizeof(real_t) );
    mass_velocity_x[1] = calloc ( (N+2),  sizeof(real_t) );

    velocity_x = calloc ( (N+2), sizeof(real_t) );
    acceleration_x = calloc ( (N+2), sizeof(real_t) );

    // Data initialization
    for ( int_t x=1; x<=N; x++ )
    {
        PN(x) = 1e-3;
        PNU(x) = 0.0;

        real_t c = x-N/2;
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
domain_save ( int_t iteration )
{
    int_t index = iteration / snapshot_frequency;
    char filename[256];
    memset ( filename, 0, 256*sizeof(char) );
    sprintf ( filename, "data/%.5ld.bin", index );

    // TODO 6 MPI I/O

    FILE *out = fopen ( filename, "wb" );
    if ( ! out ) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        exit(1);
    }
    fwrite( &mass[0][1], sizeof(real_t), N, out );
    fclose ( out );
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

