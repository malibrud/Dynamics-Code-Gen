// Validation harness: runs the generated sparse dynamics (build/dynamics.c)
// and the dense 6x6 reference (dense.h) on random states and compares qdd.
// Exit code 0 = PASS, 1 = FAIL. See dense.h for the reference implementation.

#include <stdio.h>
#include "dense.h"

// The generated implementation under test (compiled from build/dynamics.c).
void fwd_dyn_crba( double *q, double *qd, double *tau, double *fext, double *out_qdd );

int main( void ) {
    const int NTRIALS = 1000;
    const double TOL = 1e-9;
    double maxErr = 0;
    int worstTrial = -1;

    for ( int t = 0 ; t < NTRIALS ; t++ ) {
        double q[ NBOD ], qd[ NBOD ], tau[ NBOD ], fext[ NBOD * 6 ];

        // First trials pin down special cases; the rest are fully random.
        for ( int i = 0 ; i < NBOD ; i++ ) {
            q[i]   = ( t == 0 )          ? 0 : 3.0 * frand();
            qd[i]  = ( t <= 1 )          ? 0 : 2.0 * frand();
            tau[i] = ( t <= 2 )          ? 0 : 5.0 * frand();
        }
        for ( int k = 0 ; k < NBOD * 6 ; k++ )
            fext[k] = ( t <= 3 ) ? 0 : 5.0 * frand();

        double qddGen[ NBOD ], qddRef[ NBOD ];
        fwd_dyn_crba( q, qd, tau, fext, qddGen );
        fd_dense( q, qd, tau, fext, qddRef );

        for ( int i = 0 ; i < NBOD ; i++ ) {
            double err = fabs( qddGen[i] - qddRef[i] );
            if ( err > maxErr ) { maxErr = err; worstTrial = t; }
        }
    }

    printf( "validate: %d trials, max |qdd_generated - qdd_dense| = %.3e (worst at trial %d)\n",
            NTRIALS, maxErr, worstTrial );
    if ( maxErr < TOL ) {
        printf( "validate: PASS\n" );
        return 0;
    }
    printf( "validate: FAIL (tolerance %.1e)\n", TOL );
    return 1;
}
