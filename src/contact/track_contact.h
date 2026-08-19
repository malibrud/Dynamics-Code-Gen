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

    // Half-width of the road-wheel face weight's fade, as a fraction of
    // linkLen. Adjacent shoes' weights are complementary over it, so this sets
    // how far a wheel's load is shared across a joint and nothing else: at 0.25
    // the middle half of every face carries full weight and the outer quarters
    // hand over. It cannot be zero -- that is the hard in/out test the weight
    // exists to avoid.
    double faceFade;

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
    p->faceFade   = 0.25;
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

// ---------------------------------------------------------------------------
// ROAD WHEEL AGAINST THE SHOE'S INNER RUNNING FACE.
//
// The other primitive, and the one that applies to every roller NOT marked
// 'end'. A road wheel does not roll on the pins -- it rolls on the flat inner
// surface of the shoe, between them. Modelling it as pin-against-circle gets
// two things wrong, both of which showed up in practice:
//
//   - At rest the shoe faces cut R - sqrt(R^2 - ds^2/4) INTO the wheel: 9.4 mm
//     on m113 (ds 0.1504, R 0.305). The wheel is embedded in the surface it is
//     supposed to be rolling on, and the error grows as the square of pitch.
//   - A wheel resting on points can fall BETWEEN two of them once the pitch
//     reaches 2R. A wheel resting on a face cannot, at any pitch. modelLoad
//     used to warn about this for every roller; that warning was diagnosing
//     the primitive, not the model, and now covers 'end' rollers only.
//
// Worked entirely in the SHOE's frame, where the face is flat and axis-aligned:
// the running face is the plane z = runZ (= padZ + thickness) spanning
// x in [0, linkLen], and the rollers are on its +z side because the band wraps
// the outside.
//
//   c   roller axle centre, shoe coordinates
//   uc  velocity of the roller's AXLE relative to the shoe, shoe coordinates
//   us  SLIP at the contact foot: the wheel's material point there against the
//       shoe's, shoe coordinates
//
// Those are two different vectors and the law needs both. The normal damper is
// the rate conjugate of the penetration, which is set by the axle closing on
// the face (uc); the wheel's spin does not press it in. Friction is the slip of
// the surfaces in contact (us), which is what a rolling wheel drives to zero
// while its axle is moving quite fast.
//
// THE WEIGHT, and why a plain in-range test will not do. Shoes tile the band
// end to end, so a hard test on x in [0, linkLen] hands the load from one shoe
// to the next discontinuously -- the total stays put but each body's share
// jumps, which is exactly what makes an explicit integrator chatter. Nor is a
// symmetric fade enough: shoe i+1's face BEGINS where shoe i's ends, so fading
// both out near a joint would leave a support dip at every pin.
//
// What works is a partition of unity. Each face's weight is allowed to reach
// fade metres PAST both of its own ends, and the two smoothsteps are
// complementary, so for a locally straight band
//
//     w_i(x) + w_{i+1}(x - linkLen) == 1
//
// identically across the joint: at x = linkLen both are exactly 1/2. Support
// transfers smoothly and is never doubled or dropped.
//
// AND THE WEIGHT MAKES THIS A GATED POTENTIAL, U = 1/2 w(x) k delta^2, so
// -grad U has a TANGENTIAL part 1/2 k delta^2 w'(x) alongside the normal one.
// This is the same trap as the sprocket tooth law's missing g'(delta)*U0 term,
// in a new place; it is included below and src/contact/main.c checks the whole
// gradient by finite difference rather than waiting for a trajectory to notice.
// ---------------------------------------------------------------------------
// The two halves are kept apart because they act at DIFFERENT POINTS: the
// conservative part is a function of the axle's position in the shoe frame and
// so acts at the axle, while friction is a surface effect and acts at the foot,
// which is the only way it can spin the wheel.
typedef struct {
    double fCons[3];  // conservative + normal damping, at the AXLE, shoe coords
    double fFric[3];  // friction, at the FOOT, shoe coords
    double f[3];      // their sum, for sweeps and tests
    double fN;        // normal magnitude (>= 0), along -z on the shoe
    double delta;     // penetration, positive = wheel inside the face
    double w;         // partition-of-unity weight
    double foot[3];   // contact point on the face, shoe coordinates
} TrkFaceForce;

// Partition-of-unity weight and its derivative for a face spanning [0, L],
// fading over +/- fade about each end.
static double trkFaceW( double x, double L, double fade, double *dw ) {
    double a = ( x + fade ) / ( 2.0 * fade );            // rises 0 -> 1 about x = 0
    double b = ( x - L + fade ) / ( 2.0 * fade );        // rises 0 -> 1 about x = L
    double sa = a <= 0 ? 0 : a >= 1 ? 1 : a*a*( 3.0 - 2.0*a );
    double sb = b <= 0 ? 0 : b >= 1 ? 1 : b*b*( 3.0 - 2.0*b );
    double da = a <= 0 || a >= 1 ? 0 : 6.0*a*( 1.0 - a ) / ( 2.0*fade );
    double db = b <= 0 || b >= 1 ? 0 : 6.0*b*( 1.0 - b ) / ( 2.0*fade );
    if ( dw ) *dw = da * ( 1.0 - sb ) - sa * db;
    return sa * ( 1.0 - sb );
}

static int trkFaceForce( const TrackContactParams *prm, const ModelTrack *tr,
                         const ModelTrackRoller *ro, const double c[3],
                         const double uc[3], const double us[3], TrkFaceForce *out ) {
    memset( out, 0, sizeof *out );

    double runZ = tr->padZ + tr->thickness;
    double d    = c[2] - runZ;                 // axle height above the face
    // d <= 0 means the wheel has passed clean through to the outside of the
    // band. That is a broken state, not a contact, and there is no smooth force
    // that recovers from it -- which is what the settling gain ramp is for.
    if ( d <= 0 ) return 0;

    double delta = ro->radius - d;
    out->delta = delta;
    if ( delta <= 0 ) return 0;                // one-sided

    double fade = prm->faceFade * tr->linkLen;
    double dw, w = trkFaceW( c[0], tr->linkLen, fade, &dw );
    out->w = w;
    if ( w <= 0 ) return 0;

    double m    = tr->mass;
    double kn   = prm->gain * tr->kn;
    double cn   = 2.0 * prm->zeta * sqrt( kn * m );
    double kLat = prm->gain * prm->kLatFrac * tr->kn;
    double cLat = 2.0 * prm->zeta * sqrt( kLat * m );

    // delta grows as the axle descends toward the face, so d(delta)/dt = -uc[2].
    double ramp = delta < prm->deltaRamp ? delta / prm->deltaRamp : 1.0;
    double fN   = w * ( kn * delta + cn * ramp * ( -uc[2] ) );
    if ( fN < 0 ) fN = 0;
    out->fN = fN;

    out->foot[0] = c[0]; out->foot[1] = c[1]; out->foot[2] = runZ;

    // Lateral guide. Sign: c[1] is the AXLE's offset in the SHOE's frame, so a
    // roller sitting at +y means the shoe is displaced to -y and must be pulled
    // back toward +y -- the opposite sign from the pin law, where the offset is
    // measured the other way round. Its potential is 1/2 w kLat c[1]^2, so it
    // contributes to the tangential term below as well.
    double fy = w * ( kLat * c[1] + cLat * uc[1] );

    // Tangential completion of -grad( w(x) * U0 ), with
    // U0 = 1/2 kn delta^2 + 1/2 kLat c[1]^2. Without it the weight makes the
    // law non-conservative: a wheel crossing a joint enters on one shoe's
    // ramp-up and leaves on the other's ramp-down, closing a loop with net work.
    double U0 = 0.5 * kn * delta * delta + 0.5 * kLat * c[1] * c[1];
    double fx = U0 * dw;
    double fz = -fN;

    // --- friction ---------------------------------------------------------
    // The slip of the surfaces, in the face plane. Applied at the FOOT and not
    // at the axle: this is the term that spins the wheel, and on the axle it
    // would have no moment about the axle at all.
    //
    // SIGN, and it is the opposite of the pin law's. There the relative
    // velocity is the SHOE against the roller and the output is the force on
    // the shoe, so friction opposes u and the coefficient is negative. Here us
    // is the ROLLER against the shoe -- the same flip the lateral guide above
    // has -- while the output is still the force on the shoe, so friction must
    // DRAG the shoe along with the wheel and the coefficient is positive.
    //
    // Getting this backwards is not a small error, it reverses the sign of the
    // power: with force on the shoe f = sc*us and -f on the roller, the pair
    // does work at -f.us = -sc|us|^2, so a negative sc injects mu*fN*|us|
    // watts instead of dissipating them. On track_min that drove the band
    // round under its own friction and grew total energy by 5e6 J.
    double ffx = 0, ffy = 0;
    if ( prm->enableFriction && fN > 0 ) {
        double sp = sqrt( us[0]*us[0] + us[1]*us[1] );
        double sc = tr->mu * fN / sqrt( sp*sp + prm->vEps*prm->vEps );
        ffx = sc * us[0];
        ffy = sc * us[1];
    }

    out->fCons[0] = fx;  out->fCons[1] = fy;  out->fCons[2] = fz;
    out->fFric[0] = ffx; out->fFric[1] = ffy; out->fFric[2] = 0;
    for ( int k = 0 ; k < 3 ; k++ ) out->f[k] = out->fCons[k] + out->fFric[k];
    return 1;
}

// Wrench bookkeeping for one road wheel against every shoe of its band. Split
// out of trackContactFext because the geometry is the shoe's, not the roller's:
// everything below is built in the shoe frame and pushed back to base.
static void trkFaceRoller( const TrackContactParams *prm, const ModelTrack *tr,
                           const ModelTrackRoller *ro, double Xa[][6][6], double v[][6],
                           double *fext ) {
    int s = ro->body;

    double Es[3][3], cw[3], wr[3], vr[3];
    for ( int a = 0 ; a < 3 ; a++ )
    for ( int b = 0 ; b < 3 ; b++ ) Es[a][b] = Xa[s][a][b];
    xformOrigin( Xa[s], cw );
    for ( int k = 0 ; k < 3 ; k++ ) {
        wr[k] = Es[0][k]*v[s][0] + Es[1][k]*v[s][1] + Es[2][k]*v[s][2];  // omega, world
        vr[k] = Es[0][k]*v[s][3] + Es[1][k]*v[s][4] + Es[2][k]*v[s][5];  // axle velocity
    }

    for ( int i = 0 ; i < tr->shoes ; i++ ) {
        int j = tr->shoe0 + i;

        double Ej[3][3], rw[3], wj[3], vj[3];
        for ( int a = 0 ; a < 3 ; a++ )
        for ( int b = 0 ; b < 3 ; b++ ) Ej[a][b] = Xa[j][a][b];
        xformOrigin( Xa[j], rw );
        for ( int k = 0 ; k < 3 ; k++ ) {
            wj[k] = Ej[0][k]*v[j][0] + Ej[1][k]*v[j][1] + Ej[2][k]*v[j][2];
            vj[k] = Ej[0][k]*v[j][3] + Ej[1][k]*v[j][4] + Ej[2][k]*v[j][5];
        }

        // Axle in shoe coordinates.
        double rel[3], c[3];
        for ( int k = 0 ; k < 3 ; k++ ) rel[k] = cw[k] - rw[k];
        for ( int k = 0 ; k < 3 ; k++ )
            c[k] = Ej[k][0]*rel[0] + Ej[k][1]*rel[1] + Ej[k][2]*rel[2];

        // Axle velocity relative to the shoe's material point there.
        double ucw[3], uc[3];
        ucw[0] = vr[0] - ( vj[0] + wj[1]*rel[2] - wj[2]*rel[1] );
        ucw[1] = vr[1] - ( vj[1] + wj[2]*rel[0] - wj[0]*rel[2] );
        ucw[2] = vr[2] - ( vj[2] + wj[0]*rel[1] - wj[1]*rel[0] );
        for ( int k = 0 ; k < 3 ; k++ )
            uc[k] = Ej[k][0]*ucw[0] + Ej[k][1]*ucw[1] + Ej[k][2]*ucw[2];

        // A first pass fixes the foot, which the slip is then measured at. The
        // foot depends only on the axle position, so this is not an iteration.
        double runZ = tr->padZ + tr->thickness;
        double footL[3] = { c[0], c[1], runZ }, footW[3];
        for ( int k = 0 ; k < 3 ; k++ )
            footW[k] = rw[k] + Ej[0][k]*footL[0] + Ej[1][k]*footL[1] + Ej[2][k]*footL[2];

        double ar[3], aj[3], usw[3], us[3];
        for ( int k = 0 ; k < 3 ; k++ ) { ar[k] = footW[k] - cw[k]; aj[k] = footW[k] - rw[k]; }
        usw[0] = ( vr[0] + wr[1]*ar[2] - wr[2]*ar[1] ) - ( vj[0] + wj[1]*aj[2] - wj[2]*aj[1] );
        usw[1] = ( vr[1] + wr[2]*ar[0] - wr[0]*ar[2] ) - ( vj[1] + wj[2]*aj[0] - wj[0]*aj[2] );
        usw[2] = ( vr[2] + wr[0]*ar[1] - wr[1]*ar[0] ) - ( vj[2] + wj[0]*aj[1] - wj[1]*aj[0] );
        for ( int k = 0 ; k < 3 ; k++ )
            us[k] = Ej[k][0]*usw[0] + Ej[k][1]*usw[1] + Ej[k][2]*usw[2];

        TrkFaceForce ff;
        if ( !trkFaceForce( prm, tr, ro, c, uc, us, &ff ) ) continue;

        // Two application points, so two wrench pairs -- each computed once and
        // negated, the same discipline as the pin path.
        double Fc[3], Ff[3];
        for ( int k = 0 ; k < 3 ; k++ ) {
            Fc[k] = Ej[0][k]*ff.fCons[0] + Ej[1][k]*ff.fCons[1] + Ej[2][k]*ff.fCons[2];
            Ff[k] = Ej[0][k]*ff.fFric[0] + Ej[1][k]*ff.fFric[1] + Ej[2][k]*ff.fFric[2];
        }

        double mc[3], mf[3];
        mc[0] = cw[1]*Fc[2] - cw[2]*Fc[1];
        mc[1] = cw[2]*Fc[0] - cw[0]*Fc[2];
        mc[2] = cw[0]*Fc[1] - cw[1]*Fc[0];
        mf[0] = footW[1]*Ff[2] - footW[2]*Ff[1];
        mf[1] = footW[2]*Ff[0] - footW[0]*Ff[2];
        mf[2] = footW[0]*Ff[1] - footW[1]*Ff[0];

        for ( int k = 0 ; k < 3 ; k++ ) {
            fext[ j*6 + k     ] += mc[k] + mf[k];
            fext[ j*6 + 3 + k ] += Fc[k] + Ff[k];
            fext[ s*6 + k     ] -= mc[k] + mf[k];
            fext[ s*6 + 3 + k ] -= Fc[k] + Ff[k];
        }
    }
}

// Every body that receives a contact wrench must have been generated with
// external-force support, or its slots are silently never read: the force
// would vanish from the dynamics while still summing to zero in the fext
// array, so the Newton's-third-law check would keep passing. Fail loudly
// instead. Returns 0 if the model is usable.
// Largest explicit step this BAND is stable at, the counterpart of
// groundContactStabilityStep and for the same reason: the binding constraint is
// the friction regularization, not kn. What changed with the face primitive is
// which of the two binds.
//
// Under pin contact a wrapped roller engaged ~10 pins at once and they shared
// its load, so the per-pair normal force -- and hence the tangential slope
// mu*fN/vEps -- stayed low. The running face is a partition of unity whose
// weights sum to exactly 1, so ONE shoe's worth of contact now carries a whole
// road wheel. m113's per-pair load went up about tenfold with it, and with it
// the step limit: 1.11e4 N on a wheel gives mu*fN/vEps = 8.9e5 N.s/m against an
// 18 kg shoe, so dt < 2m/c = 4.1e-5 s -- tighter than the ground's 1.2e-4 and
// tighter than the 6.2e-5 the harness was choosing from the ground alone. The
// symptom was not a blow-up but a slow drip of ~0.05 J single-step energy
// rises on a run whose NET energy was still falling, which is exactly the
// signature a marginally unstable explicit damper leaves.
//
// fN is the static share over EVERY roller, not just the road wheels. Which
// rollers actually carry ground load is a property of the layout, not of the
// primitive: on a hull-wrapped belt the sprocket and idler sit on the lower run
// too, and on track_demo all four rollers are tangent to it. Counting road
// wheels only doubled track_demo's estimate of the per-pair load and quartered
// its step, taking its settle from 288k steps to 1.4M for no stability gain.
//
// BUT KNOW WHAT THIS COSTS. Counting every roller assumes every roller shares
// the ground load, and on m113 that is false: its sprocket and idler sit well
// above the road wheels and carry none of it, so the true denominator there is
// 5 per band and not 7. Measured on m113, the worst single-step energy rise
// goes 0.00029 J at dt 2.02e-5, 0.00161 J at 2.83e-5 and 0.05535 J at 6.22e-5
// -- about dt^5, far steeper than the dt^2..dt^3 of truncation, which puts the
// real boundary near 4e-5. This function now returns 5.65e-5 for m113, ABOVE
// that boundary; only the caller's 0.5 factor brings it back under. The old
// road-wheel denominator returned 4.04e-5 and the halving bought real margin.
//
// So the safety factor is no longer belt-and-braces here, it is load-bearing.
// The principled version counts the rollers actually near the lower run --
// modelLoad already has every roller's centre and radius in the root frame from
// the hull walk, so "bottom within a radius of the lowest bottom" is a few
// lines there -- and would give m113 5 and track_demo 4, which is what each
// model physically has. Models with no band return 0 and leave it to the ground.
static double trackContactStabilityStep( const Model *mdl, const TrackContactParams *prm ) {
    double mTotal = 0, gMag = 0;
    for ( int i = 0 ; i < mdl->nbod ; i++ ) mTotal += mdl->body[ i ].mass;
    for ( int k = 0 ; k < 3 ; k++ ) gMag += mdl->gravity[k] * mdl->gravity[k];
    gMag = sqrt( gMag );
    if ( gMag <= 0 || mTotal <= 0 ) return 0;

    int nRoller = 0;
    for ( int t = 0 ; t < mdl->ntrack ; t++ ) nRoller += mdl->track[t].nroller;
    if ( nRoller == 0 ) return 0;

    double fN = mTotal * gMag / nRoller;             // one roller's static share
    double worst = 0;
    for ( int t = 0 ; t < mdl->ntrack ; t++ ) {
        const ModelTrack *tr = &mdl->track[ t ];
        double kn = prm->gain * tr->kn;
        double cn = 2.0 * prm->zeta * sqrt( kn * tr->mass );
        double cT = prm->enableFriction ? tr->mu * fN / prm->vEps : 0;
        double cBody = cn + cT;                      // the weights sum to one shoe
        if ( cBody <= 0 ) continue;
        double dt = 2.0 * tr->mass / cBody;
        if ( worst == 0 || dt < worst ) worst = dt;
    }
    return worst;
}

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

            // Road wheels roll on the shoe's inner face, not on the pins; the
            // whole pin-against-circle path below applies to the sprocket and
            // idler only. See the note above trkFaceForce for why.
            if ( !ro->isEnd ) {
                trkFaceRoller( prm, tr, ro, Xa, v, fext );
                continue;
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

            // Road wheels: mirror trkFaceForce's elastic terms, in the shoe's
            // frame. U = w(x) * ( 1/2 kn delta^2 + 1/2 kLat cy^2 ).
            if ( !ro->isEnd ) {
                double runZ = tr->padZ + tr->thickness;
                double fade = prm->faceFade * tr->linkLen;
                for ( int i = 0 ; i < tr->shoes ; i++ ) {
                    int j = tr->shoe0 + i;
                    double Ej[3][3], rw[3], rel[3], cc[3];
                    for ( int a = 0 ; a < 3 ; a++ )
                    for ( int b = 0 ; b < 3 ; b++ ) Ej[a][b] = Xa[j][a][b];
                    xformOrigin( Xa[j], rw );
                    for ( int k = 0 ; k < 3 ; k++ ) rel[k] = c[k] - rw[k];
                    for ( int k = 0 ; k < 3 ; k++ )
                        cc[k] = Ej[k][0]*rel[0] + Ej[k][1]*rel[1] + Ej[k][2]*rel[2];

                    double dd = cc[2] - runZ;
                    if ( dd <= 0 ) continue;
                    double delta = ro->radius - dd;
                    if ( delta <= 0 ) continue;
                    double w = trkFaceW( cc[0], tr->linkLen, fade, NULL );
                    if ( w <= 0 ) continue;
                    U += w * ( 0.5 * kn * delta * delta + 0.5 * kLat * cc[1] * cc[1] );
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
