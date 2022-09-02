#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <math.h>

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

void domain_save ( int_t iteration );

// TODO 4: Implement the swap function to swap the content of two variables
void
swap ( m1, m2 )
{
}


int
main ( int argc, char **argv )
{
    OPTIONS *options = parse_args( argc, argv );
    if ( !options )
    {
        fprintf( stderr, "Argument parsing failed\n" );
        exit(1);
    }

    // TODO 1: Get N, max_iteration and snapshot_frequency from the options struct and store the values in the global fields with the same names


    // TODO 2: Allocate memory for the mass, velocity and acceleration arrays.
    // There should be space for N+2 elements of the real_t type in each of the arrays.
    // The arrays should also be initialized to be filled with zeroes.

    mass[0] = ;
    mass[1] = ;
    
    mass_velocity_x[0] = ;
    mass_velocity_x[1] = ;
    
    velocity_x = ;
    acceleration_x = ;

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

    for ( int_t iteration = 0; iteration <= max_iteration; iteration++ )
    {
        // TODO 3a: Update the edges of the domain based on the boundary conditions for PN and PNU.

        // TODO 3b: Update the acceleration over the entire domain and the borders

        // TODO 3c: Update the next PNU over the entire domain

        // TODO 3d: Update the next PN over the entire domain

        // TODO 3e: Update the U over the entire domain

        if ( iteration % snapshot_frequency == 0 )
        {
            printf (
                "Iteration %ld of %ld (%.2lf%% complete)\n",
                iteration,
                max_iteration,
                100.0 * (real_t) iteration / (real_t) max_iteration
            );

            domain_save ( iteration );
        }

        // TODO 4: Implement the swap function
        swap( &mass[0], &mass[1] );
        swap( &mass_velocity_x[0], &mass_velocity_x[1] );
    }

    // TODO 5: Free the heap-allocated memory


    exit ( EXIT_SUCCESS );
}


void
domain_save ( int_t iteration )
{
    int_t index = iteration / snapshot_frequency;
    char filename[256];
    memset ( filename, 0, 256*sizeof(char) );
    sprintf ( filename, "data/%.5ld.bin", index );

    FILE *out = fopen ( filename, "wb" );
    if ( ! out ) {
        fprintf(stderr, "Failed to open file: %s\n", filename);
        exit(1);
    }
    fwrite( &mass[0][1], sizeof(real_t), N, out );
    fclose ( out );
}

