// Benchmark: generated sparse dynamics (build/dynamics.c) vs the dense 6x6
// reference (src/validate/dense.h), same inputs, wall-clock timed with
// QueryPerformanceCounter. Build with /O2 -- timing unoptimized code says
// nothing useful (build.bat does this). Both take the mechanism from the
// same model file.
//
// Each implementation cycles through a pool of pre-generated random states so
// the compiler can't specialize on constant inputs, and a running checksum of
// the outputs is printed so the calls can't be dead-code eliminated. Timing
// runs in batches until a minimum elapsed time, so the per-call figure is
// stable regardless of how fast the machine is.

#include <stdio.h>
#include <windows.h>
#include "../validate/dense.h"
#include "../contact/track_contact.h"
#include "../contact/ground_contact.h"

// The generated implementation (compiled from build/dynamics.c).
void fwd_dyn_crba( double *q, double *qd, double *tau, double *fext, double *act, double *out_qdd );
extern const int fwd_dyn_crba_nq;
extern const int fwd_dyn_crba_nv;
extern const int fwd_dyn_crba_nact;

#define NSTATES 256 // power of two so the index wrap is a mask

static Model g_mdl;
static double qs[ NSTATES ][ MODEL_MAX_Q ];
static double qds[ NSTATES ][ MODEL_MAX_V ];
static double taus[ NSTATES ][ MODEL_MAX_V ];
static double fexts[ NSTATES ][ MODEL_MAX_BODIES * 6 ];
static double acts[ NSTATES ][ MODEL_MAX_ELEM ];

typedef void (*FdFn)( const double *q, const double *qd, const double *tau,
                      const double *fext, const double *act, double *qdd );

// The generated function's signature is non-const; the inputs are never
// written, so adapting via cast is safe here.
static void genAdapter( const double *q, const double *qd, const double *tau,
                        const double *fext, const double *act, double *qdd ) {
    fwd_dyn_crba( (double*)q, (double*)qd, (double*)tau, (double*)fext, (double*)act, qdd );
}

static void denseAdapter( const double *q, const double *qd, const double *tau,
                          const double *fext, const double *act, double *qdd ) {
    fd_dense( &g_mdl, q, qd, tau, fext, act, qdd );
}

// Runs fn in batches until minSeconds has elapsed; returns ns per call.
static double benchOne( FdFn fn, double minSeconds, double *checksumPerCall ) {
    LARGE_INTEGER freq, t0, t1;
    QueryPerformanceFrequency( &freq );

    double qdd[ MODEL_MAX_V ];
    double checksum = 0;
    long long calls = 0;
    double elapsed = 0;

    QueryPerformanceCounter( &t0 );
    do {
        for ( int b = 0 ; b < 10000 ; b++ ) {
            int s = (int)( calls & ( NSTATES - 1 ) );
            fn( qs[s], qds[s], taus[s], fexts[s], acts[s], qdd );
            checksum += qdd[ 0 ] + qdd[ g_mdl.nv - 1 ];
            calls++;
        }
        QueryPerformanceCounter( &t1 );
        elapsed = (double)( t1.QuadPart - t0.QuadPart ) / (double)freq.QuadPart;
    } while ( elapsed < minSeconds );

    *checksumPerCall = checksum / (double)calls;
    return elapsed * 1e9 / (double)calls;
}

// ---------------------------------------------------------------------------
// Contact force producers (src/contact). These are not part of the generated
// dynamics -- they are caller-side, and a step that uses them costs
// denseKinematics + the producers + fwd_dyn_crba -- so they are timed and
// reported separately rather than folded into the comparison above.
//
// WHY A RANDOM STATE IS A FAIR PLACE TO MEASURE THEM, even though contact is
// nominally state-dependent: both engines are all-pairs with no culling. Every
// (shoe, roller) pair and every pad corner is visited on every call, and only
// the force arithmetic is skipped when a pair is out of range. So the cost is
// set by the DECLARED counts rather than by how many pairs happen to be
// loaded, which is also what the per-pair figures below are normalized
// against. Measured against a settled m113 the difference is well under the
// run-to-run spread.
//
// Each producer is timed as (kinematics + zero + producer) and differenced
// against (kinematics + zero), so every sample runs at a different state
// instead of repeatedly re-walking one frozen configuration.
// ---------------------------------------------------------------------------
static double cXa[ MODEL_MAX_BODIES ][6][6];   // static: 590 KB, blows the stack
static double cVb[ MODEL_MAX_BODIES ][6];
static double cFext[ MODEL_MAX_BODIES * 6 ];
static TrackContactParams  g_tprm;
static GroundContactParams g_gprm;
static double g_sink;                          // keeps the calls from being elided

typedef void (*ContactFn)( int s );

static void cBase( int s ) {
    denseKinematics( &g_mdl, qs[s], qds[s], cXa, cVb );
    memset( cFext, 0, g_mdl.nbod * 6 * sizeof *cFext );
    g_sink += cXa[0][0][0];
}
static void cTrack( int s ) {
    cBase( s );
    trackContactFext( &g_mdl, &g_tprm, cXa, cVb, cFext );
    g_sink += cFext[5];
}
static void cGround( int s ) {
    double ft[6];
    cBase( s );
    groundContactFext( &g_mdl, &g_gprm, cXa, cVb, cFext, ft );
    g_sink += ft[5];
}

// Same batching discipline as benchOne, but min-of-3 because these figures are
// produced by differencing and so carry twice the noise.
static double benchContact( ContactFn fn, double minSeconds ) {
    LARGE_INTEGER freq, t0, t1;
    QueryPerformanceFrequency( &freq );
    double best = 1e30;

    for ( int rep = 0 ; rep < 3 ; rep++ ) {
        long long calls = 0;
        double elapsed = 0;
        QueryPerformanceCounter( &t0 );
        do {
            for ( int b = 0 ; b < 200 ; b++ ) {
                fn( (int)( calls & ( NSTATES - 1 ) ) );
                calls++;
            }
            QueryPerformanceCounter( &t1 );
            elapsed = (double)( t1.QuadPart - t0.QuadPart ) / (double)freq.QuadPart;
        } while ( elapsed < minSeconds );
        double ns = elapsed * 1e9 / (double)calls;
        if ( ns < best ) best = ns;
    }
    return best;
}

static void reportContact( const char *path, double nsGen ) {
    if ( !g_mdl.ntrack && !g_mdl.ground.has ) return;

    trackContactDefaults ( &g_mdl, &g_tprm );
    groundContactDefaults( &g_mdl, &g_gprm );

    int pairs = 0;
    for ( int t = 0 ; t < g_mdl.ntrack ; t++ )
        pairs += g_mdl.track[t].shoes * g_mdl.track[t].nroller;
    int pts = g_mdl.ground.npoint;

    benchContact( cBase, 0.05 );                       // warm-up
    double nsBase = benchContact( cBase,   0.30 );
    double nsTrk  = benchContact( cTrack,  0.30 ) - nsBase;
    double nsGnd  = benchContact( cGround, 0.30 ) - nsBase;
    if ( nsTrk < 0 ) nsTrk = 0;
    if ( nsGnd < 0 ) nsGnd = 0;

    double step = nsBase + nsTrk + nsGnd + nsGen;
    printf( "bench: %s: contact-in-the-loop step budget\n", path );
    printf( "bench:   denseKinematics    %8.1f ns  %5.1f%%\n", nsBase, 100.0*nsBase/step );
    if ( pairs )
        printf( "bench:   trackContactFext   %8.1f ns  %5.1f%%   %d pairs, %.1f ns/pair\n",
                nsTrk, 100.0*nsTrk/step, pairs, nsTrk / pairs );
    if ( pts )
        printf( "bench:   groundContactFext  %8.1f ns  %5.1f%%   %d points, %.1f ns/point\n",
                nsGnd, 100.0*nsGnd/step, pts, nsGnd / pts );
    printf( "bench:   fwd_dyn_crba       %8.1f ns  %5.1f%%\n", nsGen, 100.0*nsGen/step );
    printf( "bench:   full step          %8.1f ns  (contact %.1f%%, %.0f steps/s)\n",
            step, 100.0*( nsBase + nsTrk + nsGnd )/step, 1e9/step );
    if ( g_sink == 12345.6789 ) printf( "bench: (unreachable %g)\n", g_sink );
}

int main( int argc, char **argv ) {
    const char *path = argc > 1 ? argv[ 1 ] : "../models/threelink.mech";
    char err[ 256 ];
    if ( modelLoad( path, &g_mdl, err, sizeof err ) ) {
        fprintf( stderr, "bench: %s\n", err );
        return 1;
    }
    const char *unsup = modelUnsupported( &g_mdl );
    if ( unsup ) {
        fprintf( stderr, "bench: %s\n", unsup );
        return 1;
    }
    if ( fwd_dyn_crba_nq != g_mdl.nq || fwd_dyn_crba_nv != g_mdl.nv || fwd_dyn_crba_nact != g_mdl.nact ) {
        fprintf( stderr, "bench: dynamics.c was generated for nq=%d nv=%d nact=%d but '%s' has nq=%d nv=%d nact=%d -- regenerate\n",
                 fwd_dyn_crba_nq, fwd_dyn_crba_nv, fwd_dyn_crba_nact, path, g_mdl.nq, g_mdl.nv, g_mdl.nact );
        return 1;
    }

    for ( int s = 0 ; s < NSTATES ; s++ ) {
        for ( int k = 0 ; k < g_mdl.nq ; k++ )     qs[s][k]   = 3.0 * frand();
        for ( int k = 0 ; k < g_mdl.nv ; k++ )     qds[s][k]  = 2.0 * frand();
        for ( int k = 0 ; k < g_mdl.nv ; k++ )     taus[s][k] = 5.0 * frand();
        for ( int k = 0 ; k < g_mdl.nbod * 6 ; k++ ) fexts[s][k] = 5.0 * frand();
        for ( int k = 0 ; k < g_mdl.nact ; k++ )     acts[s][k]  = 5.0 * frand();
    }

    // Warm-up (page in code, settle clocks), then measure.
    double cs;
    benchOne( genAdapter,   0.05, &cs );
    benchOne( denseAdapter, 0.05, &cs );

    double csGen, csDense;
    double nsGen   = benchOne( genAdapter,   1.0, &csGen );
    double nsDense = benchOne( denseAdapter, 1.0, &csDense );

    printf( "bench: %s (%d bodies, nq %d, nv %d)\n", path, g_mdl.nbod, g_mdl.nq, g_mdl.nv );
    printf( "bench: generated sparse  %8.1f ns/call  (checksum/call %+.6f)\n", nsGen, csGen );
    printf( "bench: dense reference   %8.1f ns/call  (checksum/call %+.6f)\n", nsDense, csDense );
    printf( "bench: speedup %.1fx\n", nsDense / nsGen );
    reportContact( path, nsGen );
    return 0;
}
