// Trajectory harness: steps the generated dynamics forward in time and checks
// the invariants that only a trajectory can check.
//
// src/validate proves the generated qdd matches the dense reference at a
// state. That says nothing about whether qd is really dq/dt, whether the
// quaternion update is right, or whether a force model is passive -- all of
// which are properties of a trajectory. This harness covers that gap.
//
// Stage 0 (this file, today) runs the integrator-only check: with a
// symplectic scheme the energy error is O(dt), so halving dt must halve the
// drift. That is the same discipline the Chrono cross-check uses to separate
// truncation error from a real discrepancy -- a bug in the integrator or in
// the qd/dq-dt convention shows up as a drift that does NOT scale.
//
// Usage: sim_<model>.exe <model>.mech [dt] [--trace <file.csv>] [--stride <n>]
//                                      [--slope <deg>] [--fric anchor|reg] [--time <s>]
//
// With no --slope / --fric, a track model on ground gets TWO runs: the flat
// settle, then a parked hold on 15 degrees (see runTrackRuns).

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "integrate.h"
#include "../contact/track_contact.h"
#include "../contact/ground_contact.h"

// The generated implementation under test.
void fwd_dyn_crba( double *q, double *qd, double *tau, double *fext, double *act, double *out_qdd );
extern const int fwd_dyn_crba_nbod;
extern const int fwd_dyn_crba_nq;
extern const int fwd_dyn_crba_nv;
extern const int fwd_dyn_crba_nact;

// Seed q from the model's design pose. Only the track bands have a nonzero
// one; everything else is zeros with an identity quaternion per floating body
// (see modelLoad). This is the first consumer in the repo to actually start a
// band at its design pose -- the Chrono harness reads q0 to PLACE its bodies
// but then simulates from a scrambled state instead.
static void seedDesignPose( const Model *mdl, double *q, double *qd ) {
    memcpy( q, mdl->q0, mdl->nq * sizeof *q );
    memset( qd, 0, mdl->nv * sizeof *qd );
}

// Runs nsteps from the design pose with no applied load, and returns the
// largest |E(t) - E(0)| seen along the way. verlet selects the scheme: 0 is
// semi-implicit Euler, 1 is the drift/kick/drift Verlet composition.
//
// tau, fext and act are all zero, so the only forces are gravity and whatever
// the model's own force elements produce. Damping is what would break the
// check, so models carrying dampers are excluded by the caller.
//
// The two schemes share every line except the step itself, deliberately: the
// caller runs BOTH and asserts each one's order separately, so any difference
// between them has to come from the composition and not from the setup.
static double energyDrift( const Model *mdl, double dt, int nsteps, int verlet,
                           double *keScale ) {
    static double q[ MODEL_MAX_Q ], qd[ MODEL_MAX_V ], qdd[ MODEL_MAX_V ];
    static double qddPrev[ MODEL_MAX_V ], qdMid[ MODEL_MAX_V ];
    static double tau[ MODEL_MAX_V ], fext[ MODEL_MAX_BODIES * 6 ], act[ MODEL_MAX_ELEM ];

    seedDesignPose( mdl, q, qd );
    memset( tau,  0, sizeof tau );
    memset( fext, 0, sizeof fext );
    memset( act,  0, sizeof act );

    // A mechanism released from rest at a symmetric pose can just sit there,
    // which would pass any energy check trivially. Give it a nudge so the
    // trajectory is actually exercised.
    for ( int k = 0 ; k < mdl->nv ; k++ ) qd[k] = 0.3;

    // Seed the midpoint velocity predictor. One extra dynamics call, amortised
    // to nothing over nsteps -- and leaving qddPrev zero instead would make the
    // very first step first order, which is not a thing to do to a gate that
    // exists to measure the order.
    if ( verlet ) {
        fwd_dyn_crba( q, qd, tau, fext, act, qddPrev );
    }

    double ke0 = 0;
    double e0 = simEnergyK( mdl, q, qd, &ke0 ), worst = 0, keMax = ke0;
    for ( int s = 0 ; s < nsteps ; s++ ) {
        if ( verlet ) {
            simDrift( mdl, q, qd, 0.5 * dt );
            for ( int k = 0 ; k < mdl->nv ; k++ )
                qdMid[k] = qd[k] + qddPrev[k] * 0.5 * dt;
            fwd_dyn_crba( q, qdMid, tau, fext, act, qdd );
            memcpy( qddPrev, qdd, mdl->nv * sizeof *qdd );
            simKick( mdl, qd, qdd, dt );
            simDrift( mdl, q, qd, 0.5 * dt );
        } else {
            fwd_dyn_crba( q, qd, tau, fext, act, qdd );
            simStep( mdl, q, qd, qdd, dt );
        }
        double ke;
        double e = fabs( simEnergyK( mdl, q, qd, &ke ) - e0 );
        if ( e > worst ) worst = e;
        if ( ke > keMax ) keMax = ke;
    }
    if ( keScale ) *keScale = keMax;
    return worst;
}

// ---------------------------------------------------------------------------
// Track models: settle the band onto its rollers and the vehicle onto the
// ground, then check the trajectory invariants that only contact-in-the-loop
// can exercise.
//
// WHAT GROUND CONTACT CHANGED HERE. Every track contact is internal, so before
// ground existed the whole assembly was in free fall and dP/dt = M*g held
// identically -- which made the momentum check exact and completely blind to
// the contact model's details. Ground reactions are EXTERNAL, so that identity
// is gone and the check had to be rebuilt around the impulse the ground
// actually applied: dP = integral( M*g + F_ground ) dt. It keeps the property
// that earned it its keep, since a dropped reaction still shows up as a
// residual of order one rather than of order dt.
//
// In exchange, ground contact unlocks a stronger check that an internal-only
// band could not offer at all: STATIC EQUILIBRIUM. Once settled, the total
// ground normal load must equal the vehicle's weight. That is "it stands up",
// asserted numerically, and it is independent of every parameter in the force
// law.
//
// WHY THE GAINS RAMP. At the design pose the band wraps only the two end
// pulleys, so the five road wheels per side poke through the lower run by
// 0.17-0.25 m. At full stiffness that is order 1e6 N per side against a
// 1.1e5 N vehicle weight -- the band would be fired off the hull. Ramping the
// stiffness in (and damping it heavily at first) lets the band relax through
// the same force law that will run afterwards, with no separate "assembly"
// code path to get wrong.
// ---------------------------------------------------------------------------

// Worst penetration of any pin into any roller, and how many pairs are loaded.
static double bandWorstPenetration( const Model *mdl, const TrackContactParams *prm,
                                    const double *q, int *nLoaded ) {
    static double Xa[ MODEL_MAX_BODIES ][6][6], vb[ MODEL_MAX_BODIES ][6];
    static double qd0[ MODEL_MAX_V ];
    memset( qd0, 0, mdl->nv * sizeof *qd0 );
    denseKinematics( mdl, q, qd0, Xa, vb );

    double worst = 0;
    *nLoaded = 0;
    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[t];
        double runZ = tr->padZ + tr->thickness;
        double fade = prm->faceFade * tr->linkLen;
        for ( int r = 0 ; r < tr->nroller ; r++ ) {
            const ModelTrackRoller *ro = &tr->roller[r];
            int s = ro->body;
            double Es[3][3], c[3];
            for ( int a = 0 ; a < 3 ; a++ ) for ( int e = 0 ; e < 3 ; e++ ) Es[a][e] = Xa[s][a][e];
            xformOrigin( Xa[s], c );

            // The metric has to match the primitive. A road wheel rolls on the
            // shoe's inner running face, so its penetration is measured to that
            // face -- and its PINS are not constrained at all, sitting a chord's
            // worth OUTSIDE the pitch circle by construction. Measuring them
            // against the circle would report a large penetration for a wheel
            // resting exactly where it should.
            if ( !ro->isEnd ) {
                for ( int i = 0 ; i < tr->shoes ; i++ ) {
                    int j = tr->shoe0 + i;
                    double Ej[3][3], rw[3], rel[3], cc[3];
                    for ( int a = 0 ; a < 3 ; a++ ) for ( int e = 0 ; e < 3 ; e++ ) Ej[a][e] = Xa[j][a][e];
                    xformOrigin( Xa[j], rw );
                    for ( int k = 0 ; k < 3 ; k++ ) rel[k] = c[k] - rw[k];
                    for ( int k = 0 ; k < 3 ; k++ )
                        cc[k] = Ej[k][0]*rel[0] + Ej[k][1]*rel[1] + Ej[k][2]*rel[2];
                    double dd = cc[2] - runZ;
                    if ( dd <= 0 ) continue;
                    if ( trkFaceW( cc[0], tr->linkLen, fade, NULL ) <= 0 ) continue;
                    double delta = ro->radius - dd;
                    if ( delta > worst ) worst = delta;
                    if ( delta > 0 ) (*nLoaded)++;
                }
                continue;
            }

            for ( int i = 0 ; i < tr->shoes ; i++ ) {
                int j = tr->shoe0 + i;
                double p[3], rel[3], d[3];
                xformOrigin( Xa[j], p );
                for ( int k = 0 ; k < 3 ; k++ ) rel[k] = p[k] - c[k];
                for ( int k = 0 ; k < 3 ; k++ )
                    d[k] = Es[k][0]*rel[0] + Es[k][1]*rel[1] + Es[k][2]*rel[2];
                double delta = ro->radius - sqrt( d[0]*d[0] + d[2]*d[2] );
                if ( delta > worst ) worst = delta;
                if ( delta > 0 ) (*nLoaded)++;
            }
        }
    }
    (void)prm;
    return worst;
}

// Worst penetration of any shoe pad into the terrain, and how many pad corners
// are loaded. Mirrors bandWorstPenetration; the loaded count is what the
// realized ground damping ratio is computed from (see groundContactDefaults).
static double groundWorstPenetration( const Model *mdl, const double *q, int *nLoaded ) {
    static double Xa[ MODEL_MAX_BODIES ][6][6], vb[ MODEL_MAX_BODIES ][6];
    static double qd0[ MODEL_MAX_V ];

    double worst = 0;
    *nLoaded = 0;
    if ( !mdl->ground.has ) return 0;

    memset( qd0, 0, mdl->nv * sizeof *qd0 );
    denseKinematics( mdl, q, qd0, Xa, vb );

    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[t];
        if ( !tr->hasPad ) continue;
        for ( int i = 0 ; i < tr->shoes ; i++ ) {
            int j = tr->shoe0 + i;
            double r[3];
            xformOrigin( Xa[j], r );
            for ( int c = 0 ; c < 4 ; c++ ) {
                double tl[3], p[3], h, hx, hy;
                gndPadCorner( tr, c, tl );
                for ( int k = 0 ; k < 3 ; k++ )
                    p[k] = r[k] + Xa[j][0][k]*tl[0] + Xa[j][1][k]*tl[1] + Xa[j][2][k]*tl[2];
                if ( !terrainAt( &mdl->ground.terrain, p[0], p[1], &h, &hx, &hy ) ) continue;
                double g = h - p[2];
                if ( g > worst ) worst = g;
                if ( g > 0 ) (*nLoaded)++;
            }
        }
    }
    return worst;
}

// Optional dt override from the command line. The momentum balance below is
// asserted to be truncation rather than a modelling error, and the only honest
// way to say that is to re-run at a smaller step and watch it shrink -- the
// same discipline the loop-free energy check uses. Zero means "derive it".
static double gDtOverride = 0;

// --slope <deg> and --fric anchor|reg select ONE run of the given kind; with
// neither, build.bat's gate runs the flat settle and then the slope hold.
static double gSlopeDeg = 0;
static int    gFric = -1;                // -1 = default, 0 = regularized, 1 = anchored
static int    gSingle = 0;
static double gTime = 0;                 // --time <s>: run length; 0 = the gate's 4 s

// ---------------------------------------------------------------------------
// Optional per-step trace of the settling run.
//
// Everything the checks below assert is a scalar reduction of a trajectory:
// the worst penetration, the mean settled load, the largest single-step energy
// rise. That is the right shape for a gate, and the wrong shape for looking at
// WHY a run misbehaves -- for which the history itself is wanted. This writes
// that history as CSV, and nothing else in the harness changes.
//
// OFF UNLESS --trace IS GIVEN. build.bat runs this executable on every model as
// a gate, and a gate has no business writing files or paying for them.
//
// STRIDE, NOT EVERY STEP. m113 settles for 4 s at dt ~ 6e-5, which is 64k
// steps, and it carries 140 traced bodies (2 x (7 rollers + 63 shoes)) -- a
// full-rate dump is a few hundred MB to look at four seconds of settling. The
// default stride lands about TRACE_ROWS rows whatever dt turns out to be, so
// the same command gives a comparable file at any step size, and --stride is
// there when a specific rate is wanted.
// ---------------------------------------------------------------------------
#define TRACE_ROWS 2000

static const char *gTracePath = 0;    // --trace <file.csv>; 0 = no trace
static int gTraceStride = 0;          // --stride N; 0 = derive from TRACE_ROWS

// Column names, plus a comment line of run constants. The constants are there
// so a plot can be drawn from the file ALONE: the ground load means nothing
// without the weight it is meant to match, and a band snapshot cannot draw its
// rollers without their radii.
static void traceHeader( FILE *f, const Model *mdl, const char *path, double dt, double mTotal ) {
    double gMag = 0;
    for ( int k = 0 ; k < 3 ; k++ ) gMag += mdl->gravity[k] * mdl->gravity[k];
    fprintf( f, "# model=%s dt=%.6e mass_kg=%.6g weight_N=%.6g tracks=%d\n",
             path, dt, mTotal, mTotal * sqrt( gMag ), mdl->ntrack );
    // The terrain, so a band snapshot can draw the ground it is standing on.
    if ( mdl->ground.has && mdl->ground.terrain.kind == TERRAIN_PLANE )
        fprintf( f, "# terrain plane n=%.6g;%.6g;%.6g d=%.6g\n",
                 mdl->ground.terrain.n[0], mdl->ground.terrain.n[1],
                 mdl->ground.terrain.n[2], mdl->ground.terrain.d );
    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[t];
        fprintf( f, "# track %s shoes=%d linkLen=%.6g padZ=%.6g thickness=%.6g radii=",
                 tr->name, tr->shoes, tr->linkLen, tr->padZ, tr->thickness );
        for ( int r = 0 ; r < tr->nroller ; r++ )
            fprintf( f, "%s%.6g%s", r ? ";" : "", tr->roller[r].radius,
                     tr->roller[r].isEnd ? ( tr->roller[r].isDrive ? "d" : "e" ) : "" );
        fprintf( f, "\n" );
    }
    fprintf( f, "t,gain,hull_x,hull_y,hull_z,hull_roll_deg,hull_pitch_deg,hull_yaw_deg,"
                "fgnd_x,fgnd_y,fgnd_z,pen_band,n_band,pen_gnd,n_gnd" );
    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[t];
        for ( int r = 0 ; r < tr->nroller ; r++ )
            for ( int a = 0 ; a < 3 ; a++ )
                fprintf( f, ",%s_roller%d_%c", tr->name, r, "xyz"[a] );
        for ( int i = 0 ; i < tr->shoes ; i++ )
            for ( int a = 0 ; a < 3 ; a++ )
                fprintf( f, ",%s_shoe%d_%c", tr->name, i, "xyz"[a] );
    }
    fprintf( f, "\n" );
}

// One sample. Positions come from a fresh denseKinematics rather than from the
// Xa the step already built, because that one was evaluated at the midpoint
// POSITION for the force call -- half a step behind the state this row is
// labelled with.
static void traceRow( FILE *f, const Model *mdl, const TrackContactParams *prm,
                      const double *q, double t, double gain, const double *fGnd ) {
    static double Xa[ MODEL_MAX_BODIES ][6][6], vb[ MODEL_MAX_BODIES ][6];
    static double qd0[ MODEL_MAX_V ];
    memset( qd0, 0, mdl->nv * sizeof *qd0 );        // positions only: qd unused
    denseKinematics( mdl, q, qd0, Xa, vb );

    // The hull is the band's geometry reference, which is the body the rest of
    // the vehicle is rigidly hung off. Taken from the track rather than assumed
    // to be body 0.
    int hull = mdl->track[0].root;
    double c[3];
    xformOrigin( Xa[hull], c );

    // Xa carries the rotation as root -> body (that is how every contact metric
    // above uses it), so the body -> root matrix these angles are stated about
    // is its transpose. ZYX: R = Rz(yaw) Ry(pitch) Rx(roll).
    double R[3][3];
    for ( int a = 0 ; a < 3 ; a++ ) for ( int b = 0 ; b < 3 ; b++ ) R[a][b] = Xa[hull][b][a];
    double s2 = R[2][0] < -1.0 ? -1.0 : ( R[2][0] > 1.0 ? 1.0 : R[2][0] );
    double roll  = atan2( R[2][1], R[2][2] ) * 180.0 / M_PI;
    double pitch = asin( -s2 )              * 180.0 / M_PI;
    double yaw   = atan2( R[1][0], R[0][0] ) * 180.0 / M_PI;

    int nBand, nGnd;
    double penBand = bandWorstPenetration( mdl, prm, q, &nBand );
    double penGnd  = groundWorstPenetration( mdl, q, &nGnd );

    fprintf( f, "%.6f,%.4f,%.6f,%.6f,%.6f,%.4f,%.4f,%.4f,%.6g,%.6g,%.6g,%.6g,%d,%.6g,%d",
             t, gain, c[0], c[1], c[2], roll, pitch, yaw,
             fGnd[3], fGnd[4], fGnd[5], penBand, nBand, penGnd, nGnd );
    for ( int tr_i = 0 ; tr_i < mdl->ntrack ; tr_i++ ) {
        const ModelTrack *tr = &mdl->track[ tr_i ];
        for ( int r = 0 ; r < tr->nroller ; r++ ) {
            double p[3];
            xformOrigin( Xa[ tr->roller[r].body ], p );
            fprintf( f, ",%.6f,%.6f,%.6f", p[0], p[1], p[2] );
        }
        for ( int i = 0 ; i < tr->shoes ; i++ ) {
            double p[3];
            xformOrigin( Xa[ tr->shoe0 + i ], p );
            fprintf( f, ",%.6f,%.6f,%.6f", p[0], p[1], p[2] );
        }
    }
    fprintf( f, "\n" );
}

// Whether any pad corner of shoe j is into the terrain, from a given set of
// transforms (groundWorstPenetration's test, without its own kinematics pass).
static int shoeOnGround( const Model *mdl, const ModelTrack *tr, double Xa[][6][6], int j ) {
    double r[3];
    xformOrigin( Xa[j], r );
    for ( int corner = 0 ; corner < 4 ; corner++ ) {
        double tl[3], p[3], h, hx, hy;
        gndPadCorner( tr, corner, tl );
        for ( int k = 0 ; k < 3 ; k++ )
            p[k] = r[k] + Xa[j][0][k]*tl[0] + Xa[j][1][k]*tl[1] + Xa[j][2][k]*tl[2];
        if ( terrainAt( &mdl->ground.terrain, p[0], p[1], &h, &hx, &hy ) && h - p[2] > 0 ) return 1;
    }
    return 0;
}

// ---------------------------------------------------------------------------
// What one track run is asked to do. The default run is the flat settle every
// track model has always had; build.bat additionally gets the slope hold below
// for any model that stands on ground.
// ---------------------------------------------------------------------------
typedef struct {
    double slopeDeg;     // terrain pitch, degrees; 0 = as authored
    int    anchored;     // ground friction: 1 = anchored (holds), 0 = regularized
    int    brake;        // lock the drive sprockets (a parked vehicle)
    int    gateCreep;    // assert the settled creep rate (the slope-hold gate)
} TrackRun;

// RIGIDLY ROTATE the whole initial state, and the terrain with it, by deg
// about a y axis through the terrain under the base origin. Pitching the
// ground under a vehicle that was authored on the flat would start it buried
// at one end and hanging at the other; rotating EVERYTHING instead keeps the
// design pose exactly as close to equilibrium as it was, and leaves gravity --
// which the generated code bakes in -- untouched. A vehicle on a 15 degree
// plane under vertical gravity is a vehicle on a 15 degree slope.
//
// Every floating joint attaches to the base (modelLoad enforces it), so the
// whole state is reached through them: each gets q <- qs*q on its attitude and
// its origin carried round the axis. Everything else is relative and rides
// along. That reasoning is CHECKED rather than trusted -- every body's frame
// is compared against the rotated original through denseKinematics.
static int tiltTrackModel( Model *m, double *q, double deg ) {
    static double X0[ MODEL_MAX_BODIES ][6][6], X1[ MODEL_MAX_BODIES ][6][6];
    static double vb[ MODEL_MAX_BODIES ][6], qd0[ MODEL_MAX_V ];
    memset( qd0, 0, m->nv * sizeof *qd0 );

    double th = deg * 3.14159265358979323846 / 180.0;
    double cs = cos( th ), sn = sin( th );
    double Rs[3][3] = { { cs, 0, sn }, { 0, 1, 0 }, { -sn, 0, cs } };
    double qs[4] = { cos( 0.5*th ), 0, sin( 0.5*th ), 0 };

    double h, hx, hy;
    if ( !terrainAt( &m->ground.terrain, 0, 0, &h, &hx, &hy ) ) {
        printf( "sim: FAIL: --slope needs a terrain to tilt\n" );
        return 1;
    }
    double c[3] = { 0, 0, h };

    denseKinematics( m, q, qd0, X0, vb );

    for ( int i = 0 ; i < m->nbod ; i++ ) {
        const ModelBody *b = &m->body[ i ];
        if ( b->jtype != JOINT_FLOATING ) continue;
        double *qa = q + m->qidx[ i ], *pos = qa + 4;
        double w = qa[0], x = qa[1], y = qa[2], z = qa[3];
        // Hamilton product qs (x) q: rotate the body-to-base attitude.
        qa[0] = qs[0]*w - qs[1]*x - qs[2]*y - qs[3]*z;
        qa[1] = qs[0]*x + qs[1]*w + qs[2]*z - qs[3]*y;
        qa[2] = qs[0]*y - qs[1]*z + qs[2]*w + qs[3]*x;
        qa[3] = qs[0]*z + qs[1]*y - qs[2]*x + qs[3]*w;
        double o[3], r[3];
        for ( int k = 0 ; k < 3 ; k++ ) o[k] = b->tree[k] + pos[k] - c[k];
        for ( int k = 0 ; k < 3 ; k++ ) r[k] = Rs[k][0]*o[0] + Rs[k][1]*o[1] + Rs[k][2]*o[2];
        for ( int k = 0 ; k < 3 ; k++ ) pos[k] = c[k] + r[k] - b->tree[k];
    }

    // Plane n.p = d, carried round the same axis: n' = Rs n, and a point p on
    // it lands at c + Rs(p - c), so d' = n'.c + d - n.c.
    ModelTerrain *tn = &m->ground.terrain;
    double n1[3], nc = 0, n1c = 0;
    for ( int k = 0 ; k < 3 ; k++ ) n1[k] = Rs[k][0]*tn->n[0] + Rs[k][1]*tn->n[1] + Rs[k][2]*tn->n[2];
    for ( int k = 0 ; k < 3 ; k++ ) { nc += tn->n[k] * c[k]; n1c += n1[k] * c[k]; }
    tn->d = n1c + tn->d - nc;
    for ( int k = 0 ; k < 3 ; k++ ) tn->n[k] = n1[k];

    // The check: origin' = c + Rs(origin - c), and E' = E Rs^T (E maps base
    // to body, so rotating the body by Rs post-multiplies it by Rs^T).
    denseKinematics( m, q, qd0, X1, vb );
    double worst = 0;
    for ( int i = 0 ; i < m->nbod ; i++ ) {
        double o0[3], o1[3], d[3];
        xformOrigin( X0[i], o0 );
        xformOrigin( X1[i], o1 );
        for ( int k = 0 ; k < 3 ; k++ ) d[k] = o0[k] - c[k];
        for ( int k = 0 ; k < 3 ; k++ ) {
            double want = c[k] + Rs[k][0]*d[0] + Rs[k][1]*d[1] + Rs[k][2]*d[2];
            worst = fmax( worst, fabs( o1[k] - want ) );
        }
        for ( int a = 0 ; a < 3 ; a++ )
        for ( int k = 0 ; k < 3 ; k++ ) {
            double want = X0[i][a][0]*Rs[k][0] + X0[i][a][1]*Rs[k][1] + X0[i][a][2]*Rs[k][2];
            worst = fmax( worst, fabs( X1[i][a][k] - want ) );
        }
    }
    if ( worst > 1e-9 ) {
        printf( "sim: FAIL: tilting the initial state is not a rigid rotation (%.3g off) -- a "
                "body is reached other than through a base-attached floating joint\n", worst );
        return 1;
    }
    return 0;
}

static int runTrackSim( const Model *mdl, const char *path, const TrackRun *run ) {
    static double q[ MODEL_MAX_Q ], qd[ MODEL_MAX_V ], qdd[ MODEL_MAX_V ];
    static double tau[ MODEL_MAX_V ], fext[ MODEL_MAX_BODIES * 6 ], act[ MODEL_MAX_ELEM ];
    static double Xa[ MODEL_MAX_BODIES ][6][6], vb[ MODEL_MAX_BODIES ][6];
    static GroundAnchors anchors;
    static Model tilted;

    char err[ 256 ];
    if ( trackContactCheckModel( mdl, err, sizeof err ) ) {
        printf( "sim: FAIL: %s\n", err );
        return 1;
    }
    if ( groundContactCheckModel( mdl, err, sizeof err ) ) {
        printf( "sim: FAIL: %s\n", err );
        return 1;
    }

    memcpy( q, mdl->q0, mdl->nq * sizeof *q );
    memset( qd,  0, mdl->nv * sizeof *qd );
    memset( tau, 0, mdl->nv * sizeof *tau );
    memset( act, 0, mdl->nact * sizeof *act );

    if ( run->slopeDeg != 0 ) {
        if ( !mdl->ground.has ) {
            printf( "sim: FAIL: --slope needs a 'ground' block\n" );
            return 1;
        }
        tilted = *mdl;
        if ( tiltTrackModel( &tilted, q, run->slopeDeg ) ) return 1;
        mdl = &tilted;
    }

    TrackContactParams prm;
    trackContactDefaults( mdl, &prm );
    GroundContactParams gprm;
    groundContactDefaults( mdl, &gprm );
    memset( &anchors, 0, sizeof anchors );
    if ( run->anchored ) gprm.anchors = &anchors;

    int nLoaded0, nGnd0;
    double pen0  = bandWorstPenetration( mdl, &prm, q, &nLoaded0 );
    double gpen0 = groundWorstPenetration( mdl, q, &nGnd0 );

    // The design pose must not start meaningfully inside the ground. Asserted
    // rather than quietly corrected: a model whose ride height and ground plane
    // disagree is an authoring error, and a harness that silently lifted the
    // vehicle would hide it. Note that dropping a model ONTO the ground would
    // be just as wrong -- m113's design pose already is its ride height, and
    // lowering it until the lowest pad touched would bury the road wheels,
    // which carry no contact points, half a metre under the terrain.
    //
    // The threshold is the static deflection the model is going to settle to
    // anyway, W/(n*kn). A model authored exactly at its ride height (as
    // track_demo is) starts at a penetration of zero plus rounding, so a bare
    // "> 0" would be a rounding test rather than a modelling one.
    double sinkStatic = 0;
    if ( mdl->ground.has && mdl->ground.kn > 0 ) {
        double W = 0, gMag = 0;
        for ( int i = 0 ; i < mdl->nbod ; i++ ) W += mdl->body[i].mass;
        for ( int k = 0 ; k < 3 ; k++ ) gMag += mdl->gravity[k] * mdl->gravity[k];
        sinkStatic = W * sqrt( gMag ) / ( groundLoadedEstimate( mdl ) * mdl->ground.kn );
    }
    if ( gpen0 > sinkStatic ) {
        printf( "sim: FAIL: the design pose starts %.4g m inside the ground, past the %.4g m "
                "it would settle to -- move the model or the 'ground' plane\n", gpen0, sinkStatic );
        return 1;
    }

    // Step size. Ground friction, not stiffness, is what limits an explicit
    // integrator here (see groundContactStabilityStep), and a light-shoe model
    // like track_demo lands well under the 1e-4 the band alone was happy with.
    // Derive it rather than hard-coding one that happens to suit m113.
    // 4 s, not the 0.6 s the band alone needed. Standing the vehicle up excites
    // a slow heave/pitch mode on the torsion bars -- about 1 s per cycle on
    // m113 -- and the settled load is not meaningful until that has rung down.
    // T_RAMP IS SHORT ON PURPOSE, and it used to be 0.2 s. The ramp existed to
    // relax a band jammed 0.22 m into its own rollers by the two-pulley belt
    // idealization. Now that the belt is built around the full roller hull the
    // design pose starts within 0.009 m of equilibrium, and a long ramp does
    // the opposite of what it was for: while contact is still at low gain the
    // whole vehicle free-falls, and 0.2 s of that is 1/2*9.81*0.2^2 = 0.196 m
    // -- almost exactly the depth the road wheels were then found buried in the
    // band. The gains still ease in, over a window short enough that the drop
    // is sub-millimetre, and the artificial viscosity still drains what is left.
    const double T_RAMP = 0.02, T_SETTLE = gTime > 0 ? gTime : 4.0;
    //
    // BOTH engines get a say, and which one binds depends on the model. The
    // ground's friction sets it for a light-shoe band like track_demo; the
    // band's own face friction sets it for m113, where the partition of unity
    // puts a whole road wheel's load on one shoe. Taking only the ground's
    // estimate chose 6.2e-5 s for m113 against the band's 4.1e-5 limit, and the
    // run answered with 10469 single-step energy rises of ~0.05 J while its NET
    // energy still fell -- a marginally unstable damper, not truncation.
    //
    // FRICTION IS SOLVED IMPLICITLY HERE, so both estimates below report the
    // normal spring rather than the friction damper, and the step they return
    // is one to two orders of magnitude looser than it used to be. Setting
    // prm.dt / gprm.dt is what selects that; leaving them zero would give the
    // old explicit law and the old (much tighter) estimates. The two are set
    // twice on purpose -- once here so the estimate sees the right mode, then
    // again to the chosen DT below.
    prm.dt = gprm.dt = 1.0;                  // any positive value: mode select
    double dtGnd = groundContactStabilityStep( mdl, &gprm );
    double dtBand = trackContactStabilityStep( mdl, &prm );
    double dtStab = dtGnd;
    if ( dtBand > 0 && ( dtStab == 0 || dtBand < dtStab ) ) dtStab = dtBand;

    // THE STEP IS CHOSEN BY ACCURACY, NOT STABILITY, and that is a change worth
    // being explicit about. With friction implicit, the stability estimates
    // above come back at 1e-3 to 5e-3 s and nothing about them binds -- the run
    // would survive there and be worthless. What binds is resolving the contact
    // spring's oscillation, so the cap is a fixed number of samples per contact
    // PERIOD (see trackContactAccuracyStep for the calibration).
    //
    // This replaces a flat DT = 1e-4. That constant fit m113 and was 20% too
    // loose for track_min, which the sweep only exposed once second-order
    // stepping raised every model's usable step; under semi-implicit Euler
    // 1e-4 was inside all three models' limits and the difference never showed.
    // 90 rather than the 77-84 the sweep showed as the boundary: at 90 samples
    // m113 uses 66% of its single-step energy budget, which is about the margin
    // semi-implicit Euler had at its own 1e-4. At 80 it uses 83% and the
    // measured failure is at 68, so the default would be sitting a few percent
    // from the cliff. The override (argv[2]) is there for anyone who wants the
    // last 12%.
    const double SAMPLES = 90.0;
    double dtAccG = groundContactAccuracyStep( mdl, &gprm, SAMPLES );
    double dtAccB = trackContactAccuracyStep ( mdl, &prm,  SAMPLES );
    double DT = dtAccB;
    if ( dtAccG > 0 && ( DT == 0 || dtAccG < DT ) ) DT = dtAccG;
    if ( DT <= 0 ) DT = 1e-4;                            // no band, no pads
    // Stability still gets a veto. It does not bind on any model today, but a
    // stiffer damper or a heavier shoe could put it back under the accuracy
    // limit, and finding that out from a diverged run would be the wrong way.
    if ( dtStab > 0 && 0.5 * dtStab < DT ) DT = 0.5 * dtStab;
    if ( gDtOverride > 0 ) DT = gDtOverride;
    prm.dt = gprm.dt = DT;
    int nsteps = (int)( T_SETTLE / DT );

    // PARKING BRAKE. Without one a tracked vehicle on a slope is not parked,
    // it is coasting: the sprocket is a free revolute, the road wheels roll on
    // the band, and the hull rolls down its own track whatever the ground
    // friction does. So a slope-HOLD test has to lock the drive, and the lock
    // is a stiff torsional spring-damper on each sprocket joint -- what a
    // locked brake is, mechanically, down to its own compliance.
    //
    // Stiffness is sized off the step rather than picked: the brake's period
    // on the bare sprocket is TWO contact periods, so the step chosen above
    // resolves it with room to spare and it never becomes what sets dt.
    // Critically damped on the same inertia.
    int    nBrake = 0, brakeV[ 2 * MODEL_MAX_TRACKS ], brakeQ[ 2 * MODEL_MAX_TRACKS ];
    double brakeLock[ 2 * MODEL_MAX_TRACKS ], brakeK[ 2 * MODEL_MAX_TRACKS ], brakeC[ 2 * MODEL_MAX_TRACKS ];
    if ( run->brake ) {
        for ( int t = 0 ; t < mdl->ntrack ; t++ )
        for ( int r = 0 ; r < mdl->track[t].nroller ; r++ ) {
            const ModelTrackRoller *ro = &mdl->track[t].roller[r];
            if ( !ro->isDrive || nBrake >= 2 * MODEL_MAX_TRACKS ) continue;
            const ModelBody *b = &mdl->body[ ro->body ];
            double I = b->inertia[1] + b->mass * ( b->cg[0]*b->cg[0] + b->cg[2]*b->cg[2] );
            double w = 3.14159265358979323846 / ( SAMPLES * DT );   // 2*pi / (2 periods)
            brakeV[ nBrake ]    = mdl->vidx[ ro->body ];
            brakeQ[ nBrake ]    = mdl->qidx[ ro->body ];
            brakeLock[ nBrake ] = q[ mdl->qidx[ ro->body ] ];
            brakeK[ nBrake ]    = I * w * w;
            brakeC[ nBrake ]    = 2.0 * sqrt( brakeK[ nBrake ] * I );
            nBrake++;
        }
        if ( nBrake == 0 ) {
            printf( "sim: FAIL: a parked run needs a drive sprocket to brake, and this model has none\n" );
            return 1;
        }
        printf( "sim: %s: parking brake on %d sprocket(s), %.3g N.m/rad\n", path, nBrake, brakeK[0] );
    }
    double mTotal = 0;
    for ( int i = 0 ; i < mdl->nbod ; i++ ) mTotal += mdl->body[i].mass;

    // Trace file, if asked for. Opened here because the stride is derived from
    // nsteps, which is only known once the step size has been chosen.
    FILE *trace = 0;
    int traceStride = gTraceStride, traceRows = 0;
    // The ground wrench is AVERAGED over the stride rather than sampled. A
    // force alternating sign every step is exactly what an even stride
    // samples one phase of forever: m113's idler friction once held a +/-5 kN
    // period-2 cycle (see trkPairMass), and the sampled trace drew a steady
    // -5.8 kN on a vehicle at rest against a true mean of a few N. The cycle
    // is fixed and now gated below; the averaging stays, because it costs
    // nothing and a trace that can alias is a trace that can lie.
    double traceF[6] = { 0 };
    int traceN = 0;
    if ( gTracePath ) {
        trace = fopen( gTracePath, "w" );
        if ( !trace ) {
            printf( "sim: FAIL: cannot open trace file '%s'\n", gTracePath );
            return 1;
        }
        if ( traceStride <= 0 ) {
            traceStride = nsteps / TRACE_ROWS;
            if ( traceStride < 1 ) traceStride = 1;
        }
        traceHeader( trace, mdl, path, DT, mTotal );
    }
    if ( mdl->ground.has )
        printf( "sim: %s: ground: %d pads, kn %.3g N/m/pt, mu %.3g, mEff %.3g kg, "
                "implicit friction (mImpulse %.3g kg), stability dt %.2e s "
                "(ground %.2e, band %.2e), accuracy dt %.2e s at %.0f samples/period "
                "(ground %.2e, band %.2e) -> using dt %.2e s\n",
                path, mdl->ground.npoint, mdl->ground.kn, mdl->ground.mu, gprm.mEff,
                gprm.mImpulse, dtStab, dtGnd, dtBand,
                dtAccG > 0 && ( dtAccB == 0 || dtAccG < dtAccB ) ? dtAccG : dtAccB,
                SAMPLES, dtAccG, dtAccB, DT );

    double ePrev = 0, eFirst = 0, eLast = 0, keScale = 0;
    int    energyRises = 0;
    double worstRise = 0;

    // Momentum is only conserved once the artificial viscosity is gone: that
    // term multiplies the whole of qd, rigid-body modes included, so while it
    // is active it is an external force and dP/dt = M*g does not hold. Latch a
    // reference state at the moment it switches off and measure from there.
    double Pref[3] = { 0, 0, 0 }, tRef = 0, pRef[3] = { 0, 0, 0 };
    int    haveRef = 0;

    // Ground impulse, and its value latched at the same instant as Pref. The
    // ground wrench comes straight out of groundContactFext rather than being
    // recomputed, so the balance is measured against the very numbers the
    // dynamics was handed.
    double gImp[3] = { 0, 0, 0 }, gImpRef[3] = { 0, 0, 0 };
    double gAvg[3] = { 0, 0, 0 };
    int    nAvg = 0, avgFrom = nsteps - nsteps / 10;
    // The same window projected onto the Nyquist frequency, (-1)^s: zero for
    // anything the step resolves, and exactly A for a force alternating +/-A
    // every step. See the period-2 check after the loop.
    double gAlt[3] = { 0, 0, 0 };

    // PATCH SLIP: whether the friction is holding, measured where friction
    // acts. The shoes on the ground when the last two seconds begin are the
    // patch; their mean position is averaged over each of those seconds, and
    // the slip is the tangential part of the difference.
    //
    // NOT THE HULL. The hull sits on torsion bars, a braked sprocket and a
    // sagging top run, and it sways on them: on 15 degrees, anchored, it swings
    // +/-8 mm along the slope on a period of several seconds while the patch
    // under it moves 0.3 mm and does not drift. A hull rate measures that sway,
    // and on the flat with no brake it measures the vehicle ROLLING on its own
    // tracks -- the top run circulating at twice the hull speed over a lower
    // run that stays put. Neither is friction. Means, not end samples, because
    // the patch still breathes with the load.
    //
    // Sampled from the midpoint kinematics the force pass already computed,
    // so this costs one cross-check per shoe and no extra kinematics.
    static int patch[ MODEL_MAX_BODIES ];
    int    nPatch = 0, patchPicked = 0;
    double pWin[2][3] = { { 0 } };
    int    nWin[2] = { 0, 0 };
    const double T_WIN = 1.0;

    for ( int s = 0 ; s < nsteps ; s++ ) {
        // Drift / kick / drift, with the forces evaluated at the midpoint
        // POSITION -- the half-drift below puts q at t + DT/2. The gain ramp
        // is sampled at t + DT/2 too, since it is a coefficient of that same
        // force. The velocity is deliberately NOT advanced to the midpoint;
        // see the long note below the drift.
        double t  = s * DT;
        double tm = t + 0.5 * DT;
        double x = tm < T_RAMP ? tm / T_RAMP : 1.0;
        double sm = x*x*( 3.0 - 2.0*x );
        prm.gain = sm;
        prm.zeta = 1.0 + ( 0.2 - 1.0 ) * sm;      // over-damped early, 0.2 running
        gprm.gain = sm;
        gprm.zeta = 1.0 + ( 0.3 - 1.0 ) * sm;

        simDrift( mdl, q, qd, 0.5 * DT );

        // NO MIDPOINT VELOCITY PREDICTOR HERE, and that is the one real
        // difference between this loop and the one in energyDrift. It is worth
        // the paragraph, because the obvious "improvement" is a large
        // regression and it took a sweep to see that.
        //
        // Verlet is second order only if the force is evaluated at the MIDPOINT
        // of the step, in both slots. The position slot comes free from the
        // half-drift above. The velocity slot needs qd + qdd_prev*DT/2, and on
        // the tree models that predictor is exactly what earns the order: drop
        // it there and all eight fall from ratio 4.000 to 2.000.
        //
        // Here it does the opposite, because a predicted velocity in an
        // EXPLICIT DAMPER is an energy source. A damper -c*v evaluated at
        // v + delta but applied to v gives dv = -(c*DT/m)*(v + delta): the
        // delta part is not damping at all, it is forcing, and it does not care
        // which way v points. The implicit friction law's guarantee -- an
        // impulse that can zero the slip but never reverse it -- is stated
        // about the velocity the impulse is APPLIED to, so predicting breaks it
        // outright. Near equilibrium the true slip is on the order of vEps
        // while qdd_prev*DT/2 is set by the contact stiffness, so the predicted
        // slip can arrive with the wrong sign.
        //
        // Measured on m113, worst single-step energy rise (bound ~0.013 J):
        //
        //     dt      no predictor   predictor on both   predictor on crba only
        //     5e-5      0.000513          0.006219                 --
        //     1e-4      0.002437          0.034680             0.003175
        //     2e-4      0.010160          0.270600             0.019170
        //
        // Predicting got RELATIVELY worse as dt shrank, which is injection and
        // not truncation. The third column is why the tempting split -- predict
        // for the dynamics, whose velocity dependence is conservative Coriolis,
        // and not for contact, whose velocity dependence is dissipative -- does
        // not work either: m113's own ten torsion bars and six shocks are
        // dampers too, and they live INSIDE fwd_dyn_crba.
        //
        // So the track loop keeps the second-order POSITION and pays first
        // order in the velocity slot, which still beats semi-implicit Euler by
        // 2.4x at 1e-4 and passes at 2e-4 where Euler does not. Getting the
        // full order here means treating the contact damper and the model's
        // own dampers implicitly, the way friction already is -- a separate
        // change, and now one with a measured reason to want it.
        double fGnd[6];
        denseKinematics( mdl, q, qd, Xa, vb );
        memset( fext, 0, mdl->nbod * 6 * sizeof *fext );
        trackContactFext( mdl, &prm, Xa, vb, fext );
        groundContactFext( mdl, &gprm, Xa, vb, fext, fGnd );
        for ( int b = 0 ; b < nBrake ; b++ )
            tau[ brakeV[b] ] = -brakeK[b] * ( q[ brakeQ[b] ] - brakeLock[b] ) - brakeC[b] * qd[ brakeV[b] ];
        fwd_dyn_crba( q, qd, tau, fext, act, qdd );

        // prm.dt / gprm.dt stay at the FULL DT and that is not an oversight.
        // The implicit friction divisor k/(1 + k*dt/mEff) has to use the
        // interval over which the velocity update is applied, and the kick is
        // a whole step -- the two halves either side of it are drifts, which
        // touch no velocity at all. Halving it here to "match" the DT/2 above
        // would under-damp the friction by 2x at every contact.
        simKick( mdl, qd, qdd, DT );

        // Artificial viscosity, faded out with the ramp: it drains the
        // relaxation transient without biasing the settled state. Applied
        // between the kick and the second drift so the drift carries the
        // velocity the body actually ends the step with.
        double beta = 20.0 * ( 1.0 - sm );
        if ( beta > 0 ) { double f = 1.0 - beta*DT; for ( int k = 0 ; k < mdl->nv ; k++ ) qd[k] *= f; }

        simDrift( mdl, q, qd, 0.5 * DT );

        // Accumulated before the latch below, so gImp at latch time covers
        // exactly the steps whose forces acted before tRef. The midpoint
        // sample makes this the midpoint rule for the impulse integral, which
        // is a second-order quadrature rather than the left-hand rectangle it
        // used to be -- the momentum check gets sharper for free.
        for ( int k = 0 ; k < 3 ; k++ ) gImp[k] += fGnd[ 3 + k ] * DT;
        if ( mdl->ground.has ) {
            // Windowed on tm, the time Xa was evaluated at.
            int w = tm > T_SETTLE - T_WIN ? 1 : tm > T_SETTLE - 2.0 * T_WIN ? 0 : -1;
            if ( w >= 0 && !patchPicked ) {
                for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
                    const ModelTrack *tr = &mdl->track[t];
                    if ( !tr->hasPad ) continue;
                    for ( int i = 0 ; i < tr->shoes ; i++ )
                        if ( shoeOnGround( mdl, tr, Xa, tr->shoe0 + i ) ) patch[ nPatch++ ] = tr->shoe0 + i;
                }
                patchPicked = 1;
            }
            if ( w >= 0 && nPatch > 0 ) {
                double m[3] = { 0, 0, 0 }, o[3];
                for ( int i = 0 ; i < nPatch ; i++ ) {
                    xformOrigin( Xa[ patch[i] ], o );
                    for ( int k = 0 ; k < 3 ; k++ ) m[k] += o[k];
                }
                for ( int k = 0 ; k < 3 ; k++ ) pWin[w][k] += m[k] / nPatch;
                nWin[w]++;
            }
        }
        if ( s >= avgFrom ) {
            double sg = ( s & 1 ) ? -1.0 : 1.0;
            for ( int k = 0 ; k < 3 ; k++ ) { gAvg[k] += fGnd[ 3 + k ]; gAlt[k] += sg * fGnd[ 3 + k ]; }
            nAvg++;
        }

        // The row is labelled t + DT and written here, at the bottom of the
        // step, so its state and the impulse accumulated above refer to the
        // same instant.
        if ( trace ) {
            for ( int k = 0 ; k < 6 ; k++ ) traceF[k] += fGnd[k];
            traceN++;
        }
        if ( trace && ( s % traceStride == 0 || s == nsteps - 1 ) ) {
            for ( int k = 0 ; k < 6 ; k++ ) traceF[k] /= traceN;
            traceRow( trace, mdl, &prm, q, t + DT, sm, traceF );
            memset( traceF, 0, sizeof traceF );
            traceN = 0;
            traceRows++;
        }

        // Momentum reference, latched the instant the artificial viscosity is
        // gone -- that term multiplies the whole of qd, rigid-body modes
        // included, so while it is active it is an external force and the
        // balance below does not hold. Kept OUT of the energy window: the
        // momentum check costs nothing per step and wants the longest span it
        // can get, since it is the integrated impulse that makes a dropped
        // reaction visible.
        if ( !haveRef && t > T_RAMP * 1.5 ) {
            simMomentum( mdl, q, qd, Pref );
            for ( int k = 0 ; k < 3 ; k++ ) gImpRef[k] = gImp[k];
            tRef = t + DT;              // state is post-step
            // Where the hull was when it stopped being pushed around, for the
            // informational hull-motion line after the run. Whether FRICTION
            // held is the patch-slip measurement's job, not this one's.
            for ( int k = 0 ; k < 3 ; k++ ) pRef[k] = q[ mdl->qidx[0] + 4 + k ];
            haveRef = 1;
        }

        // Once the gains are at full strength and the artificial viscosity is
        // gone, every remaining force is either conservative or dissipative,
        // so the total energy must not increase. This is the trajectory-level
        // statement of the passivity the contact tests check pointwise, and it
        // is what a sign error or an energy-injecting damper would break.
        //
        // Restricted to the settled tail rather than the whole post-ramp
        // window, and both reasons point the same way: the budget costs a
        // second denseKinematics plus simEnergyK every step, which on 153
        // bodies dominates the run; and energy injection is easiest to see once
        // the relaxation transient is over, so the quiet quarter is the
        // stronger test as well as the cheaper one.
        if ( s >= avgFrom - nsteps / 4 && t > T_RAMP * 1.5 ) {
            // Full budget: kinetic + gravity + the model's own springs + both
            // sets of contact springs. Each term is only well defined at
            // constant gain, which is why this window is well clear of the
            // ramp -- while a gain moves, that potential is explicitly
            // time-dependent and energy is not expected to be conserved.
            static double Xe[ MODEL_MAX_BODIES ][6][6], ve[ MODEL_MAX_BODIES ][6];
            denseKinematics( mdl, q, qd, Xe, ve );
            double ke;
            double e = simEnergyK( mdl, q, qd, &ke ) + trackContactPotential( mdl, &prm, Xe )
                     + groundContactPotential( mdl, &gprm, Xe );
            for ( int b = 0 ; b < nBrake ; b++ ) {
                double dq = q[ brakeQ[b] ] - brakeLock[b];
                e += 0.5 * brakeK[b] * dq * dq;
            }
            if ( ke > keScale ) keScale = ke;
            if ( ePrev != 0 && e > ePrev ) {
                double rise = e - ePrev;
                if ( rise > worstRise ) worstRise = rise;
                energyRises++;
            }
            if ( eFirst == 0 ) eFirst = e;
            eLast = e;
            ePrev = e;
        }
    }

    if ( trace ) {
        fclose( trace );
        printf( "sim: %s: wrote %d trace rows (every %d steps, %.4g s) to %s\n",
                path, traceRows, traceStride, traceStride * DT, gTracePath );
    }

    int nLoaded1, nGnd1;
    double pen1  = bandWorstPenetration( mdl, &prm, q, &nLoaded1 );
    double gpen1 = groundWorstPenetration( mdl, q, &nGnd1 );

    printf( "sim: %s: settled %d steps (%.2f s, dt=%.0e): worst pin penetration %.4f -> %.4f m, "
            "loaded pairs %d -> %d\n",
            path, nsteps, T_SETTLE, DT, pen0, pen1, nLoaded0, nLoaded1 );

    // Static equilibrium: the check ground contact makes possible and an
    // internal-only band could not offer. Averaged over the last tenth of the
    // run because a penalty contact settles into an oscillation about the
    // static load rather than onto it.
    int gndOk = 1;
    if ( mdl->ground.has && nAvg > 0 ) {
        for ( int k = 0 ; k < 3 ; k++ ) gAvg[k] /= nAvg;
        // At rest the ground carries the whole weight vector, sum F = -M*g, so
        // comparing components is exact on a slope too and not just on the
        // flat. The lateral components are printed but NOT asserted: they are
        // the right check for slope holding, which is the next stage's job, and
        // on flat ground a band with nothing but a penalty guide horn drifts
        // slowly sideways -- track_demo shows ~10% of weight in y, which is the
        // demo's placeholder geometry rather than the ground law.
        double W = -mTotal * mdl->gravity[2];
        double relW = W != 0 ? fabs( gAvg[2] - W ) / fabs( W ) : fabs( gAvg[2] );
        // zeta_eff = zeta * sqrt( n_loaded / npoint ); see groundContactDefaults.
        double zEff = mdl->ground.npoint > 0
                    ? 0.3 * sqrt( (double)nGnd1 / mdl->ground.npoint ) : 0;
        printf( "sim: %s: ground: penetration %.4f -> %.4f m, loaded pads %d -> %d of %d "
                "(realized zeta %.3f)\n",
                path, gpen0, gpen1, nGnd0, nGnd1, mdl->ground.npoint, zEff );
        printf( "sim: %s: settled ground load %.4g N vs vehicle weight %.4g N (%.2f%% off); "
                "lateral (%.3g, %.3g) N\n",
                path, gAvg[2], W, relW * 100.0, gAvg[0], gAvg[1] );
        // 5%: the vehicle is still ringing down a 1.7 m/s landing at the end of
        // the window, and this is a settling assertion, not a precision one.
        // A vehicle that never stood up, or one held up by a sign error, misses
        // by order 100%, not by 5.
        if ( relW > 0.05 ) {
            printf( "sim: FAIL: the ground is not carrying the vehicle's weight "
                    "(%.4g N against %.4g N) -- it has not settled or is not standing up\n",
                    gAvg[2], W );
            gndOk = 0;
        }
        if ( gpen1 <= 0 ) {
            printf( "sim: FAIL: no pad is touching the ground after settling\n" );
            gndOk = 0;
        }
        // PERIOD-2. An implicit damper sized against more mass than it really
        // moves overshoots, and once the overshoot passes 2x it locks into a
        // step-to-step alternation that never decays. Nothing above sees it:
        // the mean is unaffected, the shoes move ~1e-7 m, and the energy bound
        // is relative to a kinetic energy that the alternation itself keeps
        // tiny. m113's idler friction sat in one at +/-4.8 kN for as long as
        // anyone ran it (see trkPairMass). 0.1% of weight is far above a
        // resolved signal's Nyquist leakage and far below that.
        double alt = 0;
        for ( int k = 0 ; k < 3 ; k++ ) {
            double a = fabs( gAlt[k] ) / nAvg;
            if ( a > alt ) alt = a;
        }
        printf( "sim: %s: period-2 content of the settled ground force %.3g N (%.2e of weight)\n",
                path, alt, W != 0 ? alt / fabs( W ) : alt );
        if ( alt > 1e-3 * fabs( W ) ) {
            printf( "sim: FAIL: the ground force alternates sign every step (%.4g N) -- an "
                    "implicit damper is sized against more mass than it moves\n", alt );
            gndOk = 0;
        }
    }
    if ( !gndOk ) return 1;

    // Momentum. Before ground contact every force was internal and gravity
    // uniform, so dP/dt = M*g held identically. Ground reactions are external,
    // so the target now carries their impulse as well -- measured from the
    // wrench groundContactFext actually applied, not recomputed. What the check
    // still catches is unchanged: a dropped reaction leaves a residual of order
    // one, not of order dt.
    double P[3], worstMom = 0;
    double span = nsteps * DT - tRef;
    simMomentum( mdl, q, qd, P );

    // How far the hull moved over the settled window. Informational: see the
    // note on patch slip for why this is not what says whether friction held.
    if ( mdl->ground.has && span > 0 ) {
        double d2 = 0;
        for ( int k = 0 ; k < 3 ; k++ ) {
            double dk = q[ mdl->qidx[0] + 4 + k ] - pRef[k];
            d2 += dk * dk;
        }
        printf( "sim: %s: hull moved %.4g m over the %.2f s settled window%s\n",
                path, sqrt( d2 ), span, nBrake ? "" : " (sprockets free: it can roll)" );
    }
    double slip = -1;
    if ( nPatch > 0 && nWin[0] > 0 && nWin[1] > 0 ) {
        const double *nrm = mdl->ground.terrain.n;   // plane: one normal
        double d[3], dn = 0, d2 = 0;
        for ( int k = 0 ; k < 3 ; k++ ) { d[k] = pWin[1][k] / nWin[1] - pWin[0][k] / nWin[0]; dn += d[k] * nrm[k]; }
        for ( int k = 0 ; k < 3 ; k++ ) { double tk = d[k] - dn * nrm[k]; d2 += tk * tk; }
        slip = sqrt( d2 ) / T_WIN;
        printf( "sim: %s: patch slip %.3g m/s over the last %.0f s, %d shoes on the ground "
                "(%s friction%s, %.3g deg)\n",
                path, slip, 2.0 * T_WIN, nPatch, run->anchored ? "anchored" : "regularized",
                nBrake ? ", sprockets braked" : "", run->slopeDeg );
    }
    if ( run->gateCreep ) {
        // The HOLD gate. Measured on m113 at 15 degrees, braked: anchored
        // 6.2e-5 m/s at 4 s and 1.1e-5 by 8 s, still settling and not
        // drifting; regularized a steady 1.65e-2 to 1.77e-2, the vEps creep
        // it has no way not to have. 5e-4 sits 8x above the first and 33x
        // below the second.
        if ( slip < 0 || slip > 5e-4 ) {
            printf( "sim: FAIL: a parked vehicle is not holding on the slope (patch slip %.3g m/s)\n",
                    slip );
            return 1;
        }
    }
    // Normalize against the momentum actually present, not against the
    // component's own target: gravity is z-only, so dividing by |M*g_k*span|
    // would make the x and y components an absolute test with a denominator
    // of 1 kg.m/s -- which on an 11 t vehicle is no test at all.
    double Pmag = 0;
    for ( int k = 0 ; k < 3 ; k++ ) Pmag += fabs( P[k] );
    if ( Pmag < 1.0 ) Pmag = 1.0;
    // Denominator. Pmag alone was right when every force was internal and the
    // whole vehicle was in free fall, because then P really was the largest
    // thing around. Standing on the ground it is not: the vehicle ends up
    // nearly at rest, so P collapses while the impulses that have to cancel to
    // produce that rest are each ~4e5 N.s. Dividing by Pmag then measures how
    // still the vehicle happens to have ended up, not whether a force went
    // missing -- it moved by a factor of 10 across three step sizes while the
    // absolute residual barely changed. So the scale is the larger of the
    // momentum present and the impulse actually applied. One scalar for all
    // three components, not a per-component target: gravity is z-only, so a
    // per-component denominator would leave x and y divided by ~zero.
    double impScale = 0;
    for ( int k = 0 ; k < 3 ; k++ )
        impScale += fabs( mTotal * mdl->gravity[k] * span ) + fabs( gImp[k] - gImpRef[k] );
    double scale = Pmag > impScale ? Pmag : impScale;
    double worstAbs = 0;
    for ( int k = 0 ; k < 3 ; k++ ) {
        double want = mTotal * mdl->gravity[k] * span + ( gImp[k] - gImpRef[k] );
        double got  = P[k] - Pref[k];
        double a = fabs( got - want );
        if ( a > worstAbs ) worstAbs = a;
        double e = a / scale;
        if ( e > worstMom ) worstMom = e;
    }
    printf( "sim: %s: linear momentum over the %.2f s post-ramp window matches "
            "M*g*dt + ground impulse to %.3e relative (%.4g N.s against a %.4g N.s scale)\n",
            path, span, worstMom, worstAbs, scale );

    // Band interpenetration. "It must go down" was the right assertion when the
    // only thing a settle did was relax a band jammed into its own geometry.
    // With the vehicle standing on the ground it is wrong: the load path runs
    // road wheel -> pin -> pad -> terrain, so a band that starts TANGENT to its
    // rollers (track_demo, by construction) has to compress into them to carry
    // anything at all, and 0 -> 0.0015 m is the mechanism working rather than
    // failing. What is still worth asserting is both halves of the original
    // intent, stated separately:
    //
    //   - the band must never be crushed onto its axles, loaded or not;
    //   - a band that DID start jammed must have relaxed.
    //
    // The crush threshold is a quarter of the smallest roller radius: a pin
    // that far inside the pitch circle is not in contact with a wheel in any
    // meaningful sense. m113 settles to 0.0348 m against 0.0613, having started
    // at 0.2237.
    double rMin = 0;
    for ( int t = 0 ; t < mdl->ntrack ; t++ )
        for ( int r = 0 ; r < mdl->track[t].nroller ; r++ ) {
            double rr = mdl->track[t].roller[r].radius;
            if ( rMin == 0 || rr < rMin ) rMin = rr;
        }
    double crush = 0.25 * rMin;
    if ( pen1 > crush ) {
        printf( "sim: FAIL: the band is crushed into its rollers (%.4f m against a %.4f m "
                "limit, a quarter of the smallest roller radius)\n", pen1, crush );
        return 1;
    }
    if ( pen0 > crush && pen1 >= pen0 ) {
        printf( "sim: FAIL: the band started jammed and settling did not relax it "
                "(%.4f -> %.4f m)\n", pen0, pen1 );
        return 1;
    }
    // Energy: everything left after the ramp is either conservative or
    // dissipative, so the budget must shrink. Two statements, because a
    // strictly monotone assertion would be wrong: semi-implicit Euler's energy
    // error OSCILLATES within an O(dt) band rather than decreasing, so a
    // first-order integrator on a stiff penalty law ripples upward now and
    // then no matter how correct the forces are.
    //   - net change over the window must be a loss;
    //   - no single-step rise may reach a truncation-scale fraction of the
    //     kinetic energy. A genuine energy-injecting force (a sign error, or
    //     the gated potential without its radial gradient, which is how this
    //     check earned its keep) shows up orders of magnitude above this.
    double eNet   = eLast - eFirst;
    double eBound = 1e-4 * fmax( keScale, 1.0 );
    printf( "sim: %s: energy over the post-ramp window: net %+.4g J, worst single-step rise "
            "%.4g J against a %.4g J bound (%d rises, KE scale %.4g J)\n",
            path, eNet, worstRise, eBound, energyRises, keScale );
    if ( eNet > eBound ) {
        printf( "sim: FAIL: total energy grew by %.4g J after the ramp -- "
                "the contact law is injecting energy\n", eNet );
        return 1;
    }
    if ( worstRise > eBound ) {
        printf( "sim: FAIL: a single step added %.4g J (bound %.4g J) -- "
                "that is too large to be truncation\n", worstRise, eBound );
        return 1;
    }
    // Contact is a stiff penalty law and the integrator is first order, so the
    // momentum residual is truncation, not a modelling error. That is asserted
    // rather than assumed: on m113, halving dt twice from 6.22e-5 takes the
    // absolute residual 33.4 -> 4.91 -> 1.64 N.s against a scale that stays put
    // at 8.25e5, so it is converging at least first order. A dropped reaction
    // does not converge at all -- it is order one, and here that means order
    // 1e5. Hence the tolerance can stay where it was even though the quantity
    // being normalized changed underneath it; m113 passes with ~250x margin.
    // Re-run with an explicit dt argument to repeat the scaling check.
    if ( worstMom > 1e-2 ) {
        printf( "sim: FAIL: linear momentum does not follow M*g*t + ground impulse "
                "(relative error %.3e) -- a contact reaction is being dropped\n", worstMom );
        return 1;
    }
    printf( "sim: PASS\n" );
    return 0;
}

// The track gate. A single run when the command line asked for one; otherwise
// the flat settle, then -- for a model that stands on ground -- the slope hold.
static int runTrackRuns( const Model *mdl, const char *path ) {
    const int FRIC_DEFAULT = 0;              // regularized until anchored is gated
    if ( gSingle ) {
        TrackRun one = { gSlopeDeg, gFric < 0 ? FRIC_DEFAULT : gFric, gSlopeDeg != 0, 0 };
        return runTrackSim( mdl, path, &one );
    }
    TrackRun flat = { 0, FRIC_DEFAULT, 0, 0 };
    if ( runTrackSim( mdl, path, &flat ) ) return 1;
    if ( !mdl->ground.has ) return 0;

    // THE SLOPE HOLD: parked, sprockets braked, on 15 degrees, which is well
    // inside every model's mu 0.8 (tan 15 = 0.27) so friction is asked only to
    // hold, never to slide. The trace, if any, was the flat run's.
    gTracePath = 0;
    char label[ 512 ];
    snprintf( label, sizeof label, "%s [15 deg, parked]", path );
    TrackRun hold = { 15.0, 1, 1, 1 };
    return runTrackSim( mdl, label, &hold );
}

int main( int argc, char **argv ) {
    const char *path = argc > 1 ? argv[ 1 ] : "../models/fourbar.mech";
    // argv[2] stayed a bare dt so the scaling re-runs documented above keep
    // working; the trace options are flags after it.
    for ( int i = 2 ; i < argc ; i++ ) {
        if ( !strcmp( argv[i], "--trace" ) && i + 1 < argc ) {
            gTracePath = argv[ ++i ];
        } else if ( !strcmp( argv[i], "--stride" ) && i + 1 < argc ) {
            gTraceStride = atoi( argv[ ++i ] );
        } else if ( !strcmp( argv[i], "--time" ) && i + 1 < argc ) {
            gTime = atof( argv[ ++i ] );
        } else if ( !strcmp( argv[i], "--slope" ) && i + 1 < argc ) {
            gSlopeDeg = atof( argv[ ++i ] );
            gSingle = 1;
        } else if ( !strcmp( argv[i], "--fric" ) && i + 1 < argc ) {
            const char *f = argv[ ++i ];
            if ( !strcmp( f, "anchor" ) ) gFric = 1;
            else if ( !strcmp( f, "reg" ) ) gFric = 0;
            else { fprintf( stderr, "sim: --fric takes 'anchor' or 'reg'\n" ); return 1; }
            gSingle = 1;
        } else if ( argv[i][0] != '-' ) {
            gDtOverride = atof( argv[i] );            // track models only
        } else {
            fprintf( stderr, "sim: usage: %s <model.mech> [dt] [--trace <file.csv>] "
                             "[--stride <n>] [--slope <deg>] [--fric anchor|reg] [--time <s>]\n", argv[0] );
            return 1;
        }
    }
    Model mdl;
    char err[ 256 ];
    if ( modelLoad( path, &mdl, err, sizeof err ) ) {
        fprintf( stderr, "sim: %s\n", err );
        return 1;
    }
    const char *unsup = modelUnsupported( &mdl );
    if ( unsup ) {
        fprintf( stderr, "sim: %s\n", unsup );
        return 1;
    }
    if ( fwd_dyn_crba_nbod != mdl.nbod || fwd_dyn_crba_nq != mdl.nq
      || fwd_dyn_crba_nv != mdl.nv || fwd_dyn_crba_nact != mdl.nact ) {
        fprintf( stderr, "sim: dynamics was generated for nbod=%d nq=%d nv=%d nact=%d but '%s' has\n"
                         "nbod=%d nq=%d nv=%d nact=%d -- regenerate\n",
                 fwd_dyn_crba_nbod, fwd_dyn_crba_nq, fwd_dyn_crba_nv, fwd_dyn_crba_nact,
                 path, mdl.nbod, mdl.nq, mdl.nv, mdl.nact );
        return 1;
    }

    // Track models get the contact-in-the-loop run instead of the bare
    // integrator check: they have both damping and closure loops, so the
    // energy-scaling test below does not apply to them, and settling the band
    // is the thing actually worth asserting.
    if ( mdl.ntrack ) return runTrackRuns( &mdl, path );
    if ( gTracePath )
        printf( "sim: %s: --trace covers the contact-in-the-loop track run only; this model "
                "has no band, so nothing was written\n", path );

    // Two classes of model cannot gate an energy check, for reasons that are
    // properties of the model rather than of the integrator:
    //
    //  - Damping dissipates, so energy is not conserved at all.
    //  - Loop models close their constraints with Baumgarte stabilization,
    //    which is a feedback term on constraint error: it injects and removes
    //    energy by construction and is not symplectic. Its drift is real and
    //    is NOT first order in dt, so gating on it would be gating on the
    //    stabilization scheme -- the same trap the validator avoids by
    //    sampling near the closure manifold. Loop models are exercised by the
    //    other checks instead.
    //
    // Enough tree models remain to make this a real gate, including satellite,
    // whose floating bodies are the only thing that exercises integrateQ's
    // quaternion path -- the one piece of this file that is not q += qd*dt.
    for ( int e = 0 ; e < mdl.nelem ; e++ ) {
        if ( mdl.elem[e].c != 0.0 ) {
            printf( "sim: %s: model has damping (element '%s', c = %g) -- energy check skipped\n",
                    path, mdl.elem[e].name, mdl.elem[e].c );
            return 0;
        }
    }
    if ( mdl.nloop ) {
        printf( "sim: %s: %d closure loop(s) -- Baumgarte stabilization is not symplectic,\n"
                "sim: energy check skipped (see the note in src/sim/main.c)\n", path, mdl.nloop );
        return 0;
    }

    // Same simulated duration at both step sizes, so the comparison is drift
    // per unit time and not drift per step.
    //
    // BOTH SCHEMES ARE RUN, and each is gated on its own order. That is what
    // makes this an order measurement rather than an assertion: Euler is the
    // control. If the Verlet composition were quietly still first order -- a
    // half-drift that is not self-adjoint would do it, which is exactly what
    // the exponential-map attitude update in integrate.h exists to prevent --
    // then the two would report the SAME ratio, and simply retargeting one
    // band from 2 to 4 could not tell that apart from the real thing. Running
    // both can.
    const double T = 0.5, DT = 1e-4;
    struct { const char *name; int verlet; int order; double lo, hi; } scheme[2] = {
        { "euler",  0, 1, 1.5, 4.0 },   // first order:  ratio -> 2
        { "verlet", 1, 2, 2.5, 8.0 },   // second order: ratio -> 4
    };

    int failed = 0;
    for ( int i = 0 ; i < 2 ; i++ ) {
        double keHi = 0, keLo = 0;
        double dHi = energyDrift( &mdl, DT,       (int)( T / DT ),           scheme[i].verlet, &keHi );
        double dLo = energyDrift( &mdl, DT / 2.0, (int)( T / ( DT / 2.0 ) ), scheme[i].verlet, &keLo );
        double ratio = dLo > 0 ? dHi / dLo : 0;

        // The noise floor is RELATIVE to the kinetic energy, not an absolute
        // number of joules. It used to be a flat 1e-12 J, which was fine while
        // the only scheme here was first order: nothing got near it. Verlet
        // goes straight through it -- gantry's drift lands at 4e-12 J, and at
        // 5e-5 it is LARGER than at 1e-4, which is round-off accumulating over
        // twice as many steps and not a scheme that got worse. A ratio computed
        // from round-off means nothing, so it has to be recognised as such.
        // 1e-11 of the KE scale sits about three orders above the ~1e-14
        // relative that a few thousand steps of double-precision accumulation
        // can produce, and far below any genuine truncation error.
        double floor = 1e-11 * ( keHi > keLo ? keHi : keLo );

        printf( "sim: %s: %-6s energy drift over %.2f s: %.4e at dt=%.1e, %.4e at dt=%.1e, "
                "ratio %.3f (KE scale %.3e)\n",
                path, scheme[i].name, T, dHi, DT, dLo, DT / 2.0, ratio, keHi );

        // Order p means halving dt divides the drift by 2^p. The bands are wide
        // because the drift of a symplectic scheme oscillates rather than
        // accumulating, so the max over a window is not a clean power law; what
        // this rejects is a drift that does not shrink at all (ratio ~ 1, a
        // convention or integrator bug), one that blows up, or a scheme sitting
        // at an order other than the one it claims.
        if ( dHi < floor ) {
            printf( "sim: %s: %s drift is %.2e against a %.2e noise floor -- nothing to scale\n",
                    path, scheme[i].name, dHi, floor );
            continue;
        }
        if ( ratio > scheme[i].lo && ratio < scheme[i].hi ) continue;
        printf( "sim: FAIL: %s energy drift is not order %d in dt (ratio %.3f, expected ~%d, "
                "band %.1f to %.1f)\n",
                scheme[i].name, scheme[i].order, ratio, 1 << scheme[i].order,
                scheme[i].lo, scheme[i].hi );
        failed = 1;
    }
    if ( failed ) return 1;
    printf( "sim: PASS\n" );
    return 0;
}
