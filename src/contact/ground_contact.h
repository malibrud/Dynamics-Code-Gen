// Ground contact: a body's contact points against the terrain.
//
// Sibling of track_contact.h and the same kind of thing -- a caller-side force
// producer that reads a state, sums external spatial forces into fext, and
// returns. The generated dynamics is unchanged and never learns that contact
// exists. A step with both engines is
//
//     denseKinematics ( mdl, q, qd, X, v );
//     memset( fext, 0, ... );
//     trackContactFext ( mdl, tprm, X, v, fext );   // shoes vs rollers
//     groundContactFext( mdl, gprm, X, v, fext );   // shoes vs terrain
//     fwd_dyn_crba( q, qd, tau, fext, act, qdd );
//
// The penalty-not-constraints and stateless arguments are exactly the ones in
// track_contact.h and are not repeated here.
//
// DELIBERATELY INDEPENDENT OF track_contact.h. The two laws share a damping
// formula and a regularized-Coulomb expression and nothing else -- the gate,
// the polar coordinates and the tooth potential have no counterpart here -- so
// factoring out a common header would deduplicate about four lines while
// destroying the property that a mistake in one cannot be present in the other.
// That two-independent-routes discipline is the same reason dense.h reimplements
// CRBA and gen_elems.h derives TSDA forces by a different route than the
// reference does.
//
// THE ONE STRUCTURAL DIFFERENCE from track contact: the ground is not a body.
// Every track contact is internal, so its wrenches sum to zero and the harness
// checks exactly that. Ground reactions have nowhere to go, so they are
// EXTERNAL: the total is nonzero by design, momentum is no longer conserved,
// and both harnesses had to change shape to match. See the notes in
// src/contact/main.c and src/sim/main.c.
//
// THE LAW. For a contact point at world p with world material velocity u, over
// terrain h(x,y) (see terrain.h for why the vertical gap is the strain
// measure and why grad g is left unnormalized):
//
//     g  = h(p.x,p.y) - p.z                    > 0 = penetrating
//     gd = h_x*u.x + h_y*u.y - u.z             = d(g)/dt, exactly
//     s  = kn*g + c*ramp(g)*gd,  clamped >= 0
//     F  = s * ( -h_x, -h_y, 1 )               elastic part is exactly -grad U
//     U  = 1/2 * kn * g^2
//
// Three properties fall out rather than having to be arranged:
//
//   - EXACTLY CONSERVATIVE on arbitrary terrain, by the gradient identity in
//     terrain.h. No planarity assumption anywhere.
//   - EXACTLY PASSIVE. power = F.u = -s*gd = -kn*g*gd - c*ramp*gd^2. The first
//     term is -dU/dt; the second is <= 0 unconditionally, whatever the terrain
//     does. The clamp can only zero the whole thing.
//   - NO ENGAGEMENT GATE, no polar coordinates, no atan2, no branch cut. The
//     gate in track_contact.h exists because the tooth and lateral terms are
//     gated potentials; nothing here is gated, and the one-sided g > 0 test is
//     the whole activation logic.
//
// Unlike the sprocket there is no tooth term to double-count, so friction runs
// on every ground point -- it IS the tractive force, which is the point.

#pragma once
#include <math.h>
#include <string.h>
#include "../model/kinematics.h"
#include "terrain.h"

typedef struct {
    // Damping fraction, applied against mEff below rather than against a single
    // body's mass. See the sizing note in groundContactDefaults -- this is the
    // one place the ground law must NOT copy the track law.
    double zeta;

    // Global stiffness scale, for ramped settling. Multiplies kn but not mu,
    // since friction is proportional to the normal load and is already scaled
    // through it. Same contract as TrackContactParams::gain.
    double gain;

    // Effective supported mass per contact point, for the damping derivation.
    double mEff;

    // Depth over which normal damping fades in, killing the force
    // discontinuity at first touch and sticky separation.
    double gRamp;

    // Regularization speed for Coulomb friction. Also the creep rate of a
    // parked vehicle, and -- see groundContactStabilityStep -- the quantity
    // that used to set the usable step size on its own.
    double vEps;

    // The caller's step, for SEMI-IMPLICIT FRICTION. Zero means the plain
    // explicit law, which is what every fixed-state test wants. See the note
    // above gndFricScale.
    double dt;

    // Mass the friction impulse at ONE contact point actually accelerates.
    //
    // THIS IS NOT mEff ABOVE, and the two must not be confused. mEff is the
    // SUPPORTED mass, M_total/n_loaded, which sizes the normal damper so the
    // vehicle heave mode gets the declared zeta -- it is deliberately much
    // larger than any one body. What the friction impulse accelerates is the
    // local body's share: a shoe over its four pad corners. Using the
    // supported mass here would under-limit the impulse and forfeit exactly
    // the stability guarantee the implicit form exists to provide.
    double mImpulse;

    int enableFriction;
} GroundContactParams;

// How many contact points are expected to be carrying the model, as opposed to
// how many it declares. A closed band only ever touches the ground along its
// lower run; a quarter of the perimeter is the conservative end of that (see
// the sizing note in groundContactDefaults, which explains why conservative
// means FEWER points here).
static int groundLoadedEstimate( const Model *mdl ) {
    int n = mdl->ground.npoint / 4;
    return n > 0 ? n : mdl->ground.npoint;
}

static void groundContactDefaults( const Model *mdl, GroundContactParams *p ) {
    memset( p, 0, sizeof *p );
    p->zeta  = 0.3;
    p->gain  = 1.0;
    p->gRamp = 1.0e-3;
    p->vEps  = 1.0e-2;
    p->enableFriction = 1;

    // THE EFFECTIVE MASS, and why it is not the shoe mass.
    //
    // track_contact.h derives c = 2*zeta*sqrt(k*m) with m = the shoe mass, which
    // is right there: a shoe really does bounce against a roller on its own.
    // Copying that here would be wrong, because what rests on the ground is the
    // VEHICLE, not a shoe. On m113 the mode in question is 11325 kg on a few
    // hundred point springs in parallel; damping each point against 18.02 kg
    // gives that mode zeta ~ 0.06 and it rings for seconds.
    //
    // So the share is taken against the whole model: mEff = M_total/n, with n
    // the number of points expected to be CARRYING it.
    //
    // n is not npoint. Only the lower run of a closed band ever touches the
    // ground, and by arc length that is between a quarter and a half of the
    // perimeter, so groundLoadedEstimate assumes a quarter -- the conservative
    // end for both of its uses. Assuming too few points makes mEff larger and
    // the damping stronger; it also makes the per-point load in
    // groundContactStabilityStep larger and the step estimate shorter. Both
    // errors are in the safe direction, which is the whole reason for picking
    // the low end.
    //
    // Working the heave mode through with n_loaded points actually down,
    //
    //     zeta_eff = zeta * sqrt( n_loaded / n_assumed )
    //
    // so m113 (208 loaded against 126 assumed) realizes about 1.3x the declared
    // figure. src/sim/main.c prints the loaded count and the realized zeta after
    // settling, so the assumption is checkable rather than hidden.
    //
    // Individual shoes come out over-damped by this derivation, which is
    // desirable: it kills shoe chatter.
    double mTotal = 0;
    for ( int i = 0 ; i < mdl->nbod ; i++ ) mTotal += mdl->body[ i ].mass;
    int n = groundLoadedEstimate( mdl );
    p->mEff = n > 0 ? mTotal / n : mTotal;

    // And the OTHER mass, the one the friction impulse accelerates: a shoe's
    // share over its four pad corners. Smallest across the bands, since the
    // limit has to hold for the lightest shoe in the model. Falls back to mEff
    // for a padded model with no band, which cannot happen today but would
    // otherwise leave this zero and silently disable the implicit form.
    double mIm = 0;
    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        if ( !mdl->track[t].hasPad ) continue;
        double m = mdl->track[t].mass / 4.0;
        if ( mIm == 0 || m < mIm ) mIm = m;
    }
    p->mImpulse = mIm > 0 ? mIm : p->mEff;
}

// ---------------------------------------------------------------------------
// The force law for ONE contact point against the terrain, in world
// coordinates. Kept separate from the fext bookkeeping so it can be swept
// directly, the way trkPairForce is.
//
//   p  contact point, world
//   u  material velocity of that point, world
//
// Returns 0 and leaves the output zeroed when the point is clear of the ground.
// ---------------------------------------------------------------------------
typedef struct {
    double f[3];       // total force on the body, world
    double fElast[3];  // elastic part alone; this is exactly -grad U
    double g;          // vertical gap, > 0 = penetrating
    double s;          // scalar magnitude along ( -h_x, -h_y, 1 )
    double fN;         // true normal load magnitude, = s*|grad g|
    double n[3];       // unit surface normal
} GndPointForce;

// ---------------------------------------------------------------------------
// SEMI-IMPLICIT FRICTION. Returns k, the tangential stiffness magnitude to
// apply as f = -k*u_t. Written out here rather than shared with
// track_contact.h's trkFricScale, on the same reasoning as the rest of this
// file: the two-independent-routes discipline is worth more than four
// deduplicated lines.
//
// Regularized Coulomb has slope mu*fN/vEps at zero slip, which is a viscous
// damper of that coefficient, so an explicit step needs dt < 2*m/c. That makes
// the step proportional to vEps -- and vEps is also the parked vehicle's creep
// rate, so speed could only be bought with creep. Backward Euler on
// mEff*du/dt = -k*u gives u' = u/(1 + k*dt/mEff), hence
//
//     k_implicit = k / ( 1 + k*dt/mEff )
//
// whose implied velocity change is a factor strictly inside [0, 1): friction
// can bring the slip to zero but never reverse it, at any dt. dt == 0 gives
// back the explicit law exactly, the divisor keeps |f| strictly inside the cone
// mu*fN, and it is smooth in u (an explicit impulse cap min(k*u, mEff*u/dt)
// would be stable too, but with a kink the integrator can chatter on).
//
// It does not make big steps free: the effective regularization velocity
// becomes about vEps + mu*fN*dt/mEff, so a large dt still creeps. What it buys
// is that the trade degrades smoothly instead of over a stability cliff, and
// that it is local -- lightly loaded points keep their small vEps, because the
// added term is proportional to fN.
static double gndFricScale( double mufN, double sp, double vEps, double dt, double mEff ) {
    double k = mufN / sqrt( sp*sp + vEps*vEps );
    if ( dt > 0 && mEff > 0 ) k /= 1.0 + k * dt / mEff;
    return k;
}

static int gndPointForce( const GroundContactParams *prm, const ModelGround *gnd,
                          const double p[3], const double u[3], GndPointForce *out ) {
    memset( out, 0, sizeof *out );

    double h, hx, hy;
    if ( !terrainAt( &gnd->terrain, p[0], p[1], &h, &hx, &hy ) ) return 0;

    double g = h - p[2];
    out->g = g;
    if ( g <= 0 ) return 0;                      // one-sided; no gate needed

    // grad(-g), UNNORMALIZED. Normalizing here is what would break the
    // conservativeness identity (terrain.h), and skipping it removes a sqrt.
    double G[3] = { -hx, -hy, 1.0 };
    double Gn   = sqrt( 1.0 + hx*hx + hy*hy );   // |grad g|
    for ( int k = 0 ; k < 3 ; k++ ) out->n[k] = G[k] / Gn;

    double kn = prm->gain * gnd->kn;
    double cn = 2.0 * prm->zeta * sqrt( kn * prm->mEff );

    // gd is d(g)/dt exactly -- grad(g).u -- which is what makes the damper the
    // rate conjugate of the spring and hence passive by construction.
    double gd   = hx*u[0] + hy*u[1] - u[2];
    double ramp = g < prm->gRamp ? g / prm->gRamp : 1.0;

    double sElast = kn * g;
    double s      = sElast + cn * ramp * gd;
    if ( s < 0 ) s = 0;                          // no sticky separation
    out->s  = s;
    out->fN = s * Gn;

    for ( int k = 0 ; k < 3 ; k++ ) {
        out->fElast[k] = sElast * G[k];
        out->f[k]      = s * G[k];
    }

    // --- friction ---------------------------------------------------------
    // Regularized Coulomb on the full slip in the tangent plane. No tooth term
    // competes with it here, so it runs on every point: this is the tractive
    // force. The honest cost is that a parked vehicle creeps at up to vEps.
    if ( prm->enableFriction && s > 0 && gnd->mu > 0 ) {
        double un = u[0]*out->n[0] + u[1]*out->n[1] + u[2]*out->n[2];
        double ut[3];
        for ( int k = 0 ; k < 3 ; k++ ) ut[k] = u[k] - un * out->n[k];
        double sp = sqrt( ut[0]*ut[0] + ut[1]*ut[1] + ut[2]*ut[2] );
        double sc = -gndFricScale( gnd->mu * out->fN, sp, prm->vEps,
                                   prm->dt, prm->mImpulse );
        for ( int k = 0 ; k < 3 ; k++ ) out->f[k] += sc * ut[k];
    }
    return 1;
}

// The four corners of a shoe's ground pad, in shoe-local coordinates. The pad
// is centred on the link: the shoe frame's origin is its own pin and the next
// pin is at +linkLen on x, so centring keeps the corners away from the pins,
// where two adjacent shoes' corners would otherwise land on the same world
// point and double the local stiffness while leaving gaps between.
//
// Computed on demand rather than materialized: m113 would need 63*4*2 = 504
// entries, and a new 504-entry model array is exactly the kind of thing that
// overflows the dense validator's stack.
static void gndPadCorner( const ModelTrack *tr, int corner, double t[3] ) {
    double x0 = 0.5 * ( tr->linkLen - tr->padLen );
    t[0] = x0 + ( ( corner & 1 ) ? tr->padLen : 0.0 );
    t[1] = ( ( corner & 2 ) ? 0.5 : -0.5 ) * tr->padWid;
    t[2] = tr->padZ;
}

// Every body receiving a ground wrench must have been generated with external
// force support, or its slots are silently never read. Unlike the track case
// this cannot be caught downstream at all: with no reaction body there is no
// Newton's-third-law sum to go wrong, so a dropped ground force just makes the
// vehicle fall through the floor with everything else still "passing".
// Returns 0 if the model is usable.
static int groundContactCheckModel( const Model *mdl, char *err, int errSize ) {
    if ( !mdl->ground.has ) return 0;
    if ( mdl->ground.kn <= 0 ) {
        snprintf( err, errSize, "ground: kn must be positive (got %g)", mdl->ground.kn );
        return 1;
    }
    if ( mdl->ground.terrain.kind == TERRAIN_NONE ) {
        snprintf( err, errSize, "ground: no terrain declared -- the block needs a 'plane' line" );
        return 1;
    }
    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[ t ];
        if ( tr->hasPad && !tr->hasFext ) {
            snprintf( err, errSize, "track '%s': shoes need 'fext 1' to receive ground forces "
                                    "(without it the vehicle falls through the ground with every "
                                    "other check still passing)", tr->name );
            return 1;
        }
    }
    return 0;
}

// Accumulates ground wrenches into fext (6*nbod, base coordinates, moment about
// the base origin first -- the convention in gen_rnea.h). Does NOT zero fext:
// the caller owns it and is normally also running trackContactFext.
//
// If fTotal is non-NULL it receives the summed ground wrench actually applied,
// [moment about the base origin; force]. That is not a diagnostic: it is what
// the momentum balance and the static-equilibrium check in src/sim/main.c are
// measured against, and computing it here means those checks read the same
// numbers the dynamics does rather than a second opinion about them.
static void groundContactFext( const Model *mdl, const GroundContactParams *prm,
                               double Xa[][6][6], double v[][6], double *fext,
                               double *fTotal ) {
    if ( fTotal ) memset( fTotal, 0, 6 * sizeof *fTotal );
    if ( !mdl->ground.has ) return;

    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[ t ];
        if ( !tr->hasPad ) continue;

        for ( int i = 0 ; i < tr->shoes ; i++ ) {
            int j = tr->shoe0 + i;

            double r[3];
            xformOrigin( Xa[j], r );             // shoe frame origin, world

            for ( int corner = 0 ; corner < 4 ; corner++ ) {
                double tl[3];
                gndPadCorner( tr, corner, tl );

                // World position: p = r + E^T * t, with E = Xa[j][0..2][0..2]
                // mapping base -> body.
                double p[3];
                for ( int k = 0 ; k < 3 ; k++ )
                    p[k] = r[k] + Xa[j][0][k]*tl[0] + Xa[j][1][k]*tl[1] + Xa[j][2][k]*tl[2];

                // Material velocity, body frame first: the frame origin's own
                // velocity plus omega x t. (At t = 0 this collapses to the bare
                // E^T*v[3..5] that track_contact.h uses for a pin, which is the
                // same expression seen from the shoe's own pin.)
                double ub[3];
                ub[0] = v[j][3] + v[j][1]*tl[2] - v[j][2]*tl[1];
                ub[1] = v[j][4] + v[j][2]*tl[0] - v[j][0]*tl[2];
                ub[2] = v[j][5] + v[j][0]*tl[1] - v[j][1]*tl[0];

                double u[3];
                for ( int k = 0 ; k < 3 ; k++ )
                    u[k] = Xa[j][0][k]*ub[0] + Xa[j][1][k]*ub[1] + Xa[j][2][k]*ub[2];

                GndPointForce pf;
                if ( !gndPointForce( prm, &mdl->ground, p, u, &pf ) ) continue;

                double mx = p[1]*pf.f[2] - p[2]*pf.f[1];
                double my = p[2]*pf.f[0] - p[0]*pf.f[2];
                double mz = p[0]*pf.f[1] - p[1]*pf.f[0];

                fext[ j*6 + 0 ] += mx; fext[ j*6 + 1 ] += my; fext[ j*6 + 2 ] += mz;
                fext[ j*6 + 3 ] += pf.f[0]; fext[ j*6 + 4 ] += pf.f[1]; fext[ j*6 + 5 ] += pf.f[2];

                if ( fTotal ) {
                    fTotal[0] += mx; fTotal[1] += my; fTotal[2] += mz;
                    fTotal[3] += pf.f[0]; fTotal[4] += pf.f[1]; fTotal[5] += pf.f[2];
                }
            }
        }
    }
}

// Stored elastic energy in the ground springs, for the trajectory energy
// budget. Mirrors the elastic term exactly and is meaningful only at constant
// gain, since ramping the gain deliberately makes the potential
// time-dependent -- the same caveat as trackContactPotential.
static double groundContactPotential( const Model *mdl, const GroundContactParams *prm,
                                      double Xa[][6][6] ) {
    if ( !mdl->ground.has ) return 0;
    double kn = prm->gain * mdl->ground.kn;
    double U  = 0;

    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[ t ];
        if ( !tr->hasPad ) continue;

        for ( int i = 0 ; i < tr->shoes ; i++ ) {
            int j = tr->shoe0 + i;
            double r[3];
            xformOrigin( Xa[j], r );

            for ( int corner = 0 ; corner < 4 ; corner++ ) {
                double tl[3], p[3], h, hx, hy;
                gndPadCorner( tr, corner, tl );
                for ( int k = 0 ; k < 3 ; k++ )
                    p[k] = r[k] + Xa[j][0][k]*tl[0] + Xa[j][1][k]*tl[1] + Xa[j][2][k]*tl[2];
                if ( !terrainAt( &mdl->ground.terrain, p[0], p[1], &h, &hx, &hy ) ) continue;
                double g = h - p[2];
                if ( g > 0 ) U += 0.5 * kn * g * g;
            }
        }
    }
    return U;
}

// Largest explicit step this ground model is stable at, as a diagnostic.
//
// The surprise is that kn does not set it -- the FRICTION REGULARIZATION does.
// Coulomb regularized over vEps has tangential slope mu*fN/vEps, and on m113
// that is 0.8*534/1e-2 = 4.3e4 N.s/m per point, four points on an 18 kg shoe,
// so dt < 2m/c ~ 2.1e-4 s against a normal-damping limit an order of magnitude
// looser. dt = 1e-4 clears it by only 2x, which is worth printing rather than
// discovering as a blow-up: shrinking vEps to make a parked vehicle creep less
// costs step size linearly.
//
// fN is the static share over groundLoadedEstimate points -- deliberately the
// low count, so the per-point load and hence the tangential slope come out high
// and the estimate errs short.
//
// WITH prm->dt > 0 the friction damper is solved implicitly (see gndFricScale)
// and cannot destabilize at any step, so its slope is dropped from the estimate
// and what remains is the normal spring -- which the comment above already
// notes is an order of magnitude looser. That is the whole point of the change:
// the quantity that "actually sets the usable step size" stops doing so.
static double groundContactStabilityStep( const Model *mdl, const GroundContactParams *prm ) {
    if ( !mdl->ground.has || mdl->ground.npoint <= 0 ) return 0;

    double mTotal = 0, gMag = 0;
    for ( int i = 0 ; i < mdl->nbod ; i++ ) mTotal += mdl->body[ i ].mass;
    for ( int k = 0 ; k < 3 ; k++ ) gMag += mdl->gravity[k] * mdl->gravity[k];
    gMag = sqrt( gMag );

    double kn = prm->gain * mdl->ground.kn;
    double cn = 2.0 * prm->zeta * sqrt( kn * prm->mEff );
    double fN = mTotal * gMag / groundLoadedEstimate( mdl );
    double cT = ( prm->enableFriction && prm->dt <= 0 )
                    ? mdl->ground.mu * fN / prm->vEps : 0;

    double worst = 0;
    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[ t ];
        if ( !tr->hasPad ) continue;
        double cBody = 4.0 * ( cn + cT );            // four corners on one shoe
        double dt = 2.0 * tr->mass / cBody;
        // Friction implicit: the spring's own period is what an explicit step
        // has to resolve, not its damper. Four points share one shoe, so the
        // stiffness the shoe sees is 4*kn.
        if ( prm->dt > 0 && kn > 0 ) {
            double dtOsc = 2.0 * sqrt( tr->mass / ( 4.0 * kn ) );
            if ( dtOsc < dt ) dt = dtOsc;
        }
        if ( worst == 0 || dt < worst ) worst = dt;
    }
    return worst;
}

// The step accuracy wants, the ground-side twin of trackContactAccuracyStep --
// see the note there for what the sample count means and how it was calibrated.
// Four pad corners share one shoe, so the stiffness that shoe sees is 4*kn,
// exactly as in the stability estimate above.
static double groundContactAccuracyStep( const Model *mdl, const GroundContactParams *prm,
                                         double samplesPerPeriod ) {
    if ( !mdl->ground.has || mdl->ground.npoint <= 0 || samplesPerPeriod <= 0 ) return 0;
    double kn = prm->gain * mdl->ground.kn;
    if ( kn <= 0 ) return 0;

    double worst = 0;
    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[ t ];
        if ( !tr->hasPad || tr->mass <= 0 ) continue;
        double dt = 2.0 * 3.14159265358979323846
                  * sqrt( tr->mass / ( 4.0 * kn ) ) / samplesPerPeriod;
        if ( worst == 0 || dt < worst ) worst = dt;
    }
    return worst;
}

// Convenience wrapper: kinematics + ground contact from a raw state, for
// callers that do not already have X and v (the tests). Zeroes fext.
static void groundContactFromState( const Model *mdl, const GroundContactParams *prm,
                                    const double *q, const double *qd, double *fext,
                                    double *fTotal ) {
    static double Xa[ MODEL_MAX_BODIES ][6][6];   // static: 590 KB, blows the stack
    static double vb[ MODEL_MAX_BODIES ][6];
    denseKinematics( mdl, q, qd, Xa, vb );
    memset( fext, 0, mdl->nbod * 6 * sizeof *fext );
    groundContactFext( mdl, prm, Xa, vb, fext, fTotal );
}
