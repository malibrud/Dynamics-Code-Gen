// Track-band contact forces: the band's shoes against the cylinders they ride
// on (sprocket, idler, road wheels).
//
// This is a caller-side force producer. It reads a state, writes external
// spatial forces, and returns; the generated dynamics is unchanged and never
// learns that contact exists. The step is
//
//     denseKinematics( mdl, q, qd, X, v );
//     trackContactFext( mdl, prm, X, v, fext );
//     fwd_dyn_crba( q, qd, tau, fext, act, qdd );
//
// WHY PENALTY CONTACT, and not exact unilateral constraints. The reason is
// structural, not numerical. The generated code's loop machinery is baked at
// code-generation time: nc and cidx[] are frozen in the model, the constraint
// Jacobian K is emitted as straight-line code, and the solve is sized for a
// fixed A = K H^-1 K^T. A contact set that changes every step -- pins entering
// and leaving the sprocket wrap -- would have to add and remove constraint
// rows at runtime, which cannot be done without abandoning code generation
// entirely. fext, by contrast, is a dense array read at fixed offsets, so it
// is structurally invariant. models/m113.mech says the same thing from the
// model's side: contact forces "have to be supplied from outside: apply them
// through fext on the shoes and road wheels".
//
// WHY IT IS STATELESS. trackContactFext is a pure function of (q, qd): no
// contact anchors, no active set, no memory between calls. That is what lets
// src/contact/main.c exercise it at 1000 random, off-trajectory states with
// the same machinery the rest of the project uses, and it is why the friction
// law below is regularized Coulomb rather than a bristle model. See the note
// on kt.
//
// THE PRIMITIVE. One query serves every roller: a shoe's PIN against a
// cylinder of the roller's declared radius. The shoe frame's ORIGIN IS A PIN
// (model.h: "the shoe frame is at its own pin, the next pin is at +linkLen on
// x"), so walking shoes shoe0 .. shoe0+N-1 visits all N pins exactly once with
// no local-point transform, and the pin's world velocity is just E^T*v[3..5]
// with no omega x rho term. The belt path was itself built tangent to circles
// of these radii, so the design pose is very nearly the zero-force pose.

#pragma once
#include <math.h>
#include <string.h>
#include "../model/kinematics.h"

typedef struct {
    // Damping fraction applied to every penalty term, against that term's own
    // stiffness and the shoe mass. Derived rather than declared in the .mech
    // file: it is dimensionally correct and rescales automatically when kn
    // changes. 0.2 running, 1.0 (over-damped) while settling.
    double zeta;

    // Global stiffness scale, for ramped settling. Multiplies kn, the tooth
    // stiffness and the lateral stiffness -- but NOT mu, since friction is
    // proportional to f_n and is already scaled through it.
    double gain;

    // Tooth-skip threshold: the peak circumferential force one pocket can
    // hold before the pin walks to the next one. This is a real physical
    // quantity and it is NOT in the .mech data, so it lives here. See the
    // sizing note in trackContactDefaults.
    double fSkip;

    // Distance over which normal damping fades in, killing the force
    // discontinuity at first touch and sticky separation.
    double deltaRamp;

    // Radial engagement window for the tooth and lateral terms: full strength
    // for delta > -gateOut, zero for delta < -(gateOut + gateFade).
    double gateOut, gateFade;

    // Regularization speed for Coulomb friction. A loaded stationary vehicle
    // creeps at up to about this rate -- the honest cost of staying stateless.
    double vEps;

    // Lateral (guide-horn) stiffness as a fraction of kn.
    double kLatFrac;

    // Angular position of pocket centre 0 in the sprocket's own frame. Pockets
    // sit at th0 + k*2*pi/teeth.
    double toothPhase;

    int enableTooth;     // sprocket tooth engagement
    int enableFriction;  // Coulomb on the non-drive rollers
    int enableNonEnd;    // road wheels (rollers that do not define the path)
} TrackContactParams;

static void trackContactDefaults( const Model *mdl, TrackContactParams *p ) {
    (void)mdl;
    memset( p, 0, sizeof *p );
    p->zeta       = 0.2;
    p->gain       = 1.0;
    // Sizing: at 0.3 g a ~11.3 t M113 needs ~33 kN of tractive force, ~16.6 kN
    // per track, spread over the ~5 pins on the sprocket wrap -- about 3.3 kN
    // per pin. 2e4 leaves 6x margin. Raising it makes tooth skip less likely
    // and the system stiffer (K_t = fSkip*teeth/radius).
    p->fSkip      = 2.0e4;
    p->deltaRamp  = 1.0e-3;
    p->gateOut    = 0.005;
    p->gateFade   = 0.030;
    p->vEps       = 1.0e-2;
    p->kLatFrac   = 0.1;
    p->toothPhase = 0.0;
    p->enableTooth = p->enableFriction = p->enableNonEnd = 1;
}

// C1 smoothstep engagement gate on the radial coordinate: 1 well inside the
// circle, 0 well outside, with a continuous derivative at both ends. Without
// it an explicit integrator chatters as pins enter and leave the wrap.
//
// Radial-only is enough. A pin s metres past a tangent point sits at
// rho - R = s^2/(2R), so on the m113 sprocket the gate closes at
// s = sqrt(2*0.245*0.035) = 0.131 m -- just under one shoe pitch. Departing
// pins fade out over about one link, which is the right behaviour and costs no
// design-pose data.
static double trkGate( double delta, double gateOut, double gateFade ) {
    double x = ( delta + gateOut + gateFade ) / gateFade;
    if ( x <= 0 ) return 0;
    if ( x >= 1 ) return 1;
    return x * x * ( 3.0 - 2.0 * x );
}

// d(gate)/d(delta). Needed because the gated terms are a POTENTIAL, g(delta) *
// U0(theta), and a force that uses only the circumferential part of -grad U is
// not conservative: a pin entering the wrap at one phase and leaving at another
// closes a loop with net work. The radial part is g'(delta)*U0, and it is not a
// small correction -- g' reaches 1.5/gateFade = 50 /m and the tooth potential
// reaches 2*fSkip*R/teeth = 980 J, so the omitted force is up to ~5e4 N, the
// same order as contact itself. Leaving it out pumped energy into the m113
// settle at ~13 J per step, which is how it was found.
static double trkGateSlope( double delta, double gateOut, double gateFade ) {
    double x = ( delta + gateOut + gateFade ) / gateFade;
    if ( x <= 0 || x >= 1 ) return 0;
    return 6.0 * x * ( 1.0 - x ) / gateFade;
}

// ---------------------------------------------------------------------------
// The force law for ONE pin against ONE roller, entirely in the roller's own
// coordinates: the band plane is its x-z and the axle is its y (model.h
// accepts 'plane xz' only, and every pin spins about y). Kept separate from
// the wrench bookkeeping in trackContactFext so it can be swept directly --
// the gate-smoothness and tooth tests in src/contact/main.c drive this
// function, not a whole assembled band.
//
//   d  pin centre relative to the axle, roller coordinates
//   u  velocity of the pin's material point relative to the roller's, at the
//      pin, roller coordinates
//   f  resulting force ON THE SHOE, roller coordinates
//
// Returns 0 and leaves f zeroed when the pair is out of range.
// ---------------------------------------------------------------------------
typedef struct {
    double f[3];      // total force on the shoe, roller coordinates
    double fN;        // radial magnitude (>= 0)
    double fT;        // circumferential (tooth) magnitude, signed
    double delta;     // penetration, positive = pin inside the circle
    double gate;      // radial engagement gate
} TrkPairForce;

static int trkPairForce( const TrackContactParams *prm, const ModelTrack *tr,
                         const ModelTrackRoller *ro, const double d[3], const double u[3],
                         TrkPairForce *out ) {
    memset( out, 0, sizeof *out );

    double R   = ro->radius;
    double rho = sqrt( d[0]*d[0] + d[2]*d[2] );
    if ( rho < 1e-12 ) return 0;                 // pin on the axle: no direction

    double delta = R - rho;
    double g     = trkGate( delta, prm->gateOut, prm->gateFade );
    double gp    = trkGateSlope( delta, prm->gateOut, prm->gateFade );
    out->delta = delta;
    out->gate  = g;
    if ( delta <= 0 && g <= 0 ) return 0;

    double m    = tr->mass;
    double kn   = prm->gain * tr->kn;
    double kLat = prm->gain * prm->kLatFrac * tr->kn;
    double cn   = 2.0 * prm->zeta * sqrt( kn   * m );
    double cLat = 2.0 * prm->zeta * sqrt( kLat * m );

    // Unit radial and circumferential directions in the band plane.
    double nx = d[0]/rho, nz = d[2]/rho;
    double tx = -nz,      tz = nx;

    double un = u[0]*nx + u[2]*nz;               // radial closing rate (+ = separating)
    double ut = u[0]*tx + u[2]*tz;               // circumferential slip
    double uy = u[1];                            // lateral slip

    // --- radial -----------------------------------------------------------
    // One-sided: the band wraps the OUTSIDE, so a pin at rho > R has lifted
    // off. The clamp stops the damper pulling the shoe back as it separates.
    double fN = 0;
    if ( delta > 0 ) {
        double ramp = delta < prm->deltaRamp ? delta / prm->deltaRamp : 1.0;
        fN = kn * delta + cn * ramp * ( -un );
        if ( fN < 0 ) fN = 0;
    }
    out->fN = fN;

    // The gated terms are a potential U = g(delta)*U0, so the radial part of
    // -grad U is g'(delta)*U0. Accumulated here and added to the radial force
    // below; kept OUT of fN, which is the true contact normal load that
    // friction is proportional to and that the design-pose test asserts on.
    double gradU0 = 0;

    double fx = 0, fy = 0, fz = 0;

    // --- tooth engagement -------------------------------------------------
    // f = -g * [ fSkip*sin(teeth*(th - th0)) + c*ut ].
    //
    // Because 'teeth' is an INTEGER the sine is automatically periodic at the
    // pocket pitch 2*pi/teeth and automatically continuous across th = +/-pi:
    // no phase wrapping, no round(), no branch. It vanishes at the pocket
    // centre (stable) and at the tooth tip (unstable -- correct, a pin
    // balanced on a tip has no circumferential force), so the taper at half
    // pitch is built in, and |f| <= fSkip makes tooth skip fall out of the
    // law rather than being modelled separately. Two-sided by construction,
    // which is what a pin sitting between two flanks needs to carry drive AND
    // overrun.
    if ( prm->enableTooth && ro->isDrive && ro->teeth > 0 ) {
        double kTooth = prm->gain * prm->fSkip * ro->teeth / R;  // small-offset slope
        double cTooth = 2.0 * prm->zeta * sqrt( kTooth * m );
        double th     = atan2( d[2], d[0] );
        double amp    = prm->gain * prm->fSkip;
        double phase  = ro->teeth * ( th - prm->toothPhase );
        // U0 = (amp*R/teeth) * (1 - cos(phase)). In polar coordinates the
        // conservative tangential force is -(1/rho) dU0/dtheta, so the moment
        // arm is rho and NOT R -- dividing by R instead leaves a residual of
        // order delta/rho, which reaches 14% at the edge of the gate band.
        // Using rho costs nothing and keeps the peak exactly fSkip ON the
        // pitch circle, where rho == R, which is where the unit test measures.
        double fT = -g * ( amp * ( R / rho ) * sin( phase ) + cTooth * ut );
        out->fT = fT;
        fx += fT * tx; fz += fT * tz;
        gradU0 += ( amp * R / ro->teeth ) * ( 1.0 - cos( phase ) );
    }

    // --- lateral guide ----------------------------------------------------
    // The band is a floating planar chain parented to the BASE, closed only by
    // a point-on-line loop constraining two in-plane translations. Its
    // rigid-body y-translation, roll and yaw are therefore completely
    // unconstrained, and gravity is z-only: with purely in-plane contact those
    // three modes would receive no force at all and the band would walk off
    // the wheels. The guide horns are the real mechanism and d[1] is already
    // in hand. Outside the friction cone -- a horn is a stop, not friction.
    fy += -g * ( kLat * d[1] + cLat * uy );
    gradU0 += 0.5 * kLat * d[1] * d[1];

    // Radial completion of -grad( g(delta) * U0 ), which is what makes the
    // gated elastic terms conservative. Sign: delta = R - rho, so
    // d(delta)/d(rho) = -1 and -dU/d(rho) = +g'(delta)*U0 outward.
    double fRad = fN + gp * gradU0;
    fx += fRad * nx; fz += fRad * nz;

    // --- friction ---------------------------------------------------------
    // Skipped on the drive sprocket: the tooth term already IS the
    // circumferential force there, and Coulomb on top would double-count
    // drive.
    if ( prm->enableFriction && fN > 0 && !ro->isDrive ) {
        double tux = u[0] - un*nx, tuy = u[1], tuz = u[2] - un*nz;
        double sp  = sqrt( tux*tux + tuy*tuy + tuz*tuz );
        double sc  = -tr->mu * fN / sqrt( sp*sp + prm->vEps*prm->vEps );
        fx += sc*tux; fy += sc*tuy; fz += sc*tuz;
    }

    out->f[0] = fx; out->f[1] = fy; out->f[2] = fz;
    return 1;
}

// Every body that receives a contact wrench must have been generated with
// external-force support, or its slots are silently never read: the force
// would vanish from the dynamics while still summing to zero in the fext
// array, so the Newton's-third-law check would keep passing. Fail loudly
// instead. Returns 0 if the model is usable.
static int trackContactCheckModel( const Model *mdl, char *err, int errSize ) {
    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[ t ];
        if ( !tr->hasFext ) {
            snprintf( err, errSize, "track '%s': shoes need 'fext 1' to receive contact forces", tr->name );
            return 1;
        }
        for ( int r = 0 ; r < tr->nroller ; r++ ) {
            int b = tr->roller[ r ].body;
            if ( !mdl->body[ b ].hasFext ) {
                snprintf( err, errSize,
                          "track '%s': roller body %d needs 'fext 1' to receive the contact reaction "
                          "(without it the reaction is dropped and the vehicle gains momentum from nothing)",
                          tr->name, b );
                return 1;
            }
        }
    }
    return 0;
}

// Accumulates contact wrenches into fext (6*nbod, base coordinates, moment
// about the base origin first -- the convention in gen_rnea.h). Does NOT zero
// fext: the caller owns it and may be adding other loads.
static void trackContactFext( const Model *mdl, const TrackContactParams *prm,
                              double Xa[][6][6], double v[][6], double *fext ) {
    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[ t ];

        for ( int r = 0 ; r < tr->nroller ; r++ ) {
            const ModelTrackRoller *ro = &tr->roller[ r ];
            if ( !ro->isEnd && !prm->enableNonEnd ) continue;

            int s = ro->body;

            // Roller frame: E maps base -> roller coordinates, c is its axle.
            double Es[3][3], c[3];
            for ( int a = 0 ; a < 3 ; a++ )
            for ( int b = 0 ; b < 3 ; b++ ) Es[a][b] = Xa[s][a][b];
            xformOrigin( Xa[s], c );

            double ws[3], vs[3];   // roller angular velocity and axle velocity, world
            for ( int k = 0 ; k < 3 ; k++ ) {
                ws[k] = Es[0][k]*v[s][0] + Es[1][k]*v[s][1] + Es[2][k]*v[s][2];
                vs[k] = Es[0][k]*v[s][3] + Es[1][k]*v[s][4] + Es[2][k]*v[s][5];
            }

            for ( int i = 0 ; i < tr->shoes ; i++ ) {
                int j = tr->shoe0 + i;

                double p[3];
                xformOrigin( Xa[j], p );            // the pin, exactly

                double rel[3], d[3];
                for ( int k = 0 ; k < 3 ; k++ ) rel[k] = p[k] - c[k];
                for ( int k = 0 ; k < 3 ; k++ )
                    d[k] = Es[k][0]*rel[0] + Es[k][1]*rel[1] + Es[k][2]*rel[2];

                // Relative velocity of the two material points at p, then
                // rotated into roller coordinates. The pin term needs no
                // omega x rho: p IS body j's frame origin, so its material
                // point's velocity is exactly E^T * v[j][3..5].
                double uw[3], u[3];
                for ( int k = 0 ; k < 3 ; k++ )
                    uw[k] = Xa[j][0][k]*v[j][3] + Xa[j][1][k]*v[j][4] + Xa[j][2][k]*v[j][5];
                uw[0] -= vs[0] + ws[1]*rel[2] - ws[2]*rel[1];
                uw[1] -= vs[1] + ws[2]*rel[0] - ws[0]*rel[2];
                uw[2] -= vs[2] + ws[0]*rel[1] - ws[1]*rel[0];
                for ( int k = 0 ; k < 3 ; k++ )
                    u[k] = Es[k][0]*uw[0] + Es[k][1]*uw[1] + Es[k][2]*uw[2];

                TrkPairForce pf;
                if ( !trkPairForce( prm, tr, ro, d, u, &pf ) ) continue;

                // Back to base coordinates: F_world = E^T * f_roller.
                double F[3];
                for ( int k = 0 ; k < 3 ; k++ )
                    F[k] = Es[0][k]*pf.f[0] + Es[1][k]*pf.f[1] + Es[2][k]*pf.f[2];

                // --- apply ----------------------------------------------
                // [ p x F ; F ] on the shoe, exactly negated on the roller.
                // The moment is computed ONCE, from the same p and the same F,
                // and negated, so each pair cancels EXACTLY -- the total then
                // cancels only to accumulation roundoff, since the two halves
                // land in different running sums and floating-point addition
                // is not associative. Do not "simplify" this into two cross
                // products about p and p-c: that would make even the pairwise
                // cancellation inexact and would silently move a moment onto
                // the hull through the roller's joint.
                double mx = p[1]*F[2] - p[2]*F[1];
                double my = p[2]*F[0] - p[0]*F[2];
                double mz = p[0]*F[1] - p[1]*F[0];

                fext[ j*6 + 0 ] += mx; fext[ j*6 + 1 ] += my; fext[ j*6 + 2 ] += mz;
                fext[ j*6 + 3 ] += F[0]; fext[ j*6 + 4 ] += F[1]; fext[ j*6 + 5 ] += F[2];

                fext[ s*6 + 0 ] -= mx; fext[ s*6 + 1 ] -= my; fext[ s*6 + 2 ] -= mz;
                fext[ s*6 + 3 ] -= F[0]; fext[ s*6 + 4 ] -= F[1]; fext[ s*6 + 5 ] -= F[2];
            }
        }
    }
}

// Stored elastic energy in the contact springs, for the trajectory energy
// budget: without it, a pin compressing a contact spring looks like energy
// leaving the system and releasing it looks like energy appearing. Mirrors
// trkPairForce's elastic terms exactly -- radial 1/2*kn*delta^2 plus the gated
// tooth and lateral potentials -- and is meaningful only at constant gain,
// since ramping the gain deliberately makes the potential time-dependent.
static double trackContactPotential( const Model *mdl, const TrackContactParams *prm,
                                     double Xa[][6][6] ) {
    double U = 0;
    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[ t ];
        double kn   = prm->gain * tr->kn;
        double kLat = prm->gain * prm->kLatFrac * tr->kn;

        for ( int r = 0 ; r < tr->nroller ; r++ ) {
            const ModelTrackRoller *ro = &tr->roller[ r ];
            if ( !ro->isEnd && !prm->enableNonEnd ) continue;
            int s = ro->body;
            double Es[3][3], c[3];
            for ( int a = 0 ; a < 3 ; a++ )
            for ( int b = 0 ; b < 3 ; b++ ) Es[a][b] = Xa[s][a][b];
            xformOrigin( Xa[s], c );

            for ( int i = 0 ; i < tr->shoes ; i++ ) {
                int j = tr->shoe0 + i;
                double p[3], rel[3], d[3];
                xformOrigin( Xa[j], p );
                for ( int k = 0 ; k < 3 ; k++ ) rel[k] = p[k] - c[k];
                for ( int k = 0 ; k < 3 ; k++ )
                    d[k] = Es[k][0]*rel[0] + Es[k][1]*rel[1] + Es[k][2]*rel[2];

                double rho = sqrt( d[0]*d[0] + d[2]*d[2] );
                if ( rho < 1e-12 ) continue;
                double delta = ro->radius - rho;
                double g = trkGate( delta, prm->gateOut, prm->gateFade );
                if ( delta <= 0 && g <= 0 ) continue;

                if ( delta > 0 ) U += 0.5 * kn * delta * delta;

                double u0 = 0.5 * kLat * d[1] * d[1];
                if ( prm->enableTooth && ro->isDrive && ro->teeth > 0 ) {
                    double amp   = prm->gain * prm->fSkip;
                    double phase = ro->teeth * ( atan2( d[2], d[0] ) - prm->toothPhase );
                    u0 += ( amp * ro->radius / ro->teeth ) * ( 1.0 - cos( phase ) );
                }
                U += g * u0;
            }
        }
    }
    return U;
}

// Convenience wrapper: kinematics + contact from a raw state, for callers that
// do not already have X and v (the tests). Zeroes fext.
static void trackContactFromState( const Model *mdl, const TrackContactParams *prm,
                                   const double *q, const double *qd, double *fext ) {
    static double Xa[ MODEL_MAX_BODIES ][6][6];   // static: 590 KB, blows the stack
    static double vb[ MODEL_MAX_BODIES ][6];
    denseKinematics( mdl, q, qd, Xa, vb );
    memset( fext, 0, mdl->nbod * 6 * sizeof *fext );
    trackContactFext( mdl, prm, Xa, vb, fext );
}
