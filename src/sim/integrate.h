// Plain-C time stepping for the generated forward dynamics.
//
// Until now nothing in build.bat stepped a trajectory: src/validate and
// src/bench evaluate fwd_dyn_crba at random states, and the only integrator in
// the repo lived inside src/chrono/main.cpp, which needs the Chrono SDK to
// compile. The track contact engine cannot be developed or tested without a
// stepper, so integrateQ is ported here with no dependencies beyond the model.
//
// TWO SCHEMES LIVE HERE, and both are symplectic -- energy oscillates within a
// bounded band instead of growing secularly the way explicit Euler does.
//
//   simStep   semi-implicit (symplectic) Euler, FIRST order. Velocity is
//             updated first and the NEW velocity advances the position:
//
//                 qd += qdd*dt;   q = drift( q, qd, dt )
//
//             This is the scheme stepGenerated() uses in the Chrono
//             cross-check. It is kept as the CONTROL for the order gate in
//             src/sim/main.c, which asserts one scheme against the other.
//
//   drift/kick/drift   Stormer-Verlet, SECOND order, composed by the caller:
//
//                 q = drift( q, qd, dt/2 )
//                 qdd = f( q, qd + qdd_prev*dt/2 )      <- ONE dynamics call
//                 qd += qdd*dt
//                 q = drift( q, qd, dt/2 )
//
//             Why the caller composes it rather than a simStep-shaped
//             function doing it: the forces have to be evaluated in the MIDDLE
//             of the step, and between the kinematics and the dynamics the
//             contact engine runs. A callback would be the only function
//             pointer in the repo, so simDrift and simKick are exposed instead
//             and each caller writes the five lines out.
//
//             The force evaluation sits at the midpoint, which is the reason
//             for this form rather than the algebraically equivalent velocity-
//             Verlet (kick/drift/kick): every damper here, the contact normal
//             damping and the friction law are all VELOCITY dependent, and
//             kick/drift/kick wants the force at the end state whose velocity
//             is the one being solved for. At the midpoint the position is
//             exact and the velocity comes from the previous step's
//             acceleration, which costs nothing.
//
//             THE VELOCITY PREDICTOR IS THE CALLER'S CALL, and the two callers
//             in src/sim/main.c genuinely differ. It is what makes the scheme
//             second order -- drop it and the tree models fall from an energy
//             drift ratio of 4.000 to 2.000 -- but a predicted velocity fed to
//             an EXPLICIT DAMPER is an energy source, because the delta between
//             the predicted and the updated velocity enters as forcing rather
//             than as damping. The contact-in-the-loop run therefore omits it
//             and keeps only the second-order position; see the note beside
//             its drift. That asymmetry is a property of the force model, not
//             of the integrator, which is why it lives there and not here.
//
// WHY THE DRIFT HAS TO BE EXACT. Verlet is second order because it is a
// SYMMETRIC composition, and that needs the half-drift to be self-adjoint. For
// a single-DOF joint q += qd*tau already is, exactly. For a FLOATING body it
// is not automatic, and this is the subtlety integrateQ exists for: qd carries
// the body's spatial velocity in its OWN frame while q carries a base-frame
// quaternion plus a base-frame position, so the rates are
//
//     Qdot = 0.5 * Q x (0, w_body)        rdot = R(Q) * v_body
//
// Held at constant qd the first of those has a closed form -- the exponential
// map -- so the attitude drift is written EXACTLY rather than to first order
// with a renormalisation afterwards, and the position drift uses the MID-drift
// attitude (the midpoint rule, self-adjoint and locally O(tau^3)) rather than
// the pre-step one. The earlier first-order pair was fine under Euler and
// would have quietly capped the order of the composition above.
//
// Getting this wrong is invisible to src/validate (which never integrates) and
// to the acceleration cross-check (which compares at a fixed state) -- exactly
// the blind spot that once hid the E = R(q) convention bug, so the convention
// note in model.h is the contract here too.

#pragma once
#include <math.h>
#include <string.h>
#include "../model/kinematics.h"

// Hamilton product, (w,x,y,z) order as everywhere else in the repo. The body-
// frame increment POST-multiplies: Qdot = 0.5 * Q x (0, w_body) is a right
// multiplication, so the flow it generates is Q(tau) = Q0 x exp(0.5*w*tau).
static void simQuatMul( const double a[4], const double b[4], double out[4] ) {
    double w = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
    double x = a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2];
    double y = a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1];
    double z = a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0];
    out[0] = w; out[1] = x; out[2] = y; out[3] = z;
}

// dq = exp( 0.5 * w_body * t ) = ( cos(th), sin(th)*w_hat ) with th = |w|*t/2.
// The exact attitude increment for a constant body-frame angular velocity, and
// a unit quaternion by construction -- which is what retires the renormalise-
// after-truncating step the first-order update needed.
//
// Written through sin(th)/|w| rather than sin(th)*w_hat so the |w| -> 0 limit
// is finite: below the branch the two-term series is already exact to double
// precision (th < 1e-8 leaves th^4/120 far under the ulp of 1).
static void simRotIncrement( const double *wb, double t, double dq[4] ) {
    double a  = 0.5 * t;
    double s2 = wb[0]*wb[0] + wb[1]*wb[1] + wb[2]*wb[2];
    double th = sqrt( s2 ) * a;
    double sc = th < 1e-8 ? a * ( 1.0 - th*th/6.0 ) : sin( th ) / sqrt( s2 );
    dq[0] = cos( th );
    dq[1] = sc * wb[0];
    dq[2] = sc * wb[1];
    dq[3] = sc * wb[2];
}

// q_out = q advanced by qd over dt. q_out may alias q.
static void simIntegrateQ( const Model *mdl, const double *q, const double *qd,
                           double dt, double *qOut ) {
    for ( int i = 0 ; i < mdl->nbod ; i++ ) {
        const ModelBody *b = &mdl->body[ i ];
        const double *qs = q  + mdl->qidx[ i ];
        const double *vs = qd + mdl->vidx[ i ];
        double *out = qOut + mdl->qidx[ i ];

        if ( b->jtype != JOINT_FLOATING ) {
            out[0] = qs[0] + vs[0] * dt;
            continue;
        }

        // Attitude: the EXACT flow, Q(dt) = Q0 x exp( 0.5 * w_body * dt ).
        // Applied as two half-increments rather than one whole one, which is
        // free (they share an axis, so they commute and compose exactly) and
        // hands the position update below the mid-drift attitude it wants.
        double qn[4] = { qs[0], qs[1], qs[2], qs[3] };
        double n = 1.0 / sqrt( qn[0]*qn[0] + qn[1]*qn[1] + qn[2]*qn[2] + qn[3]*qn[3] );
        for ( int k = 0 ; k < 4 ; k++ ) qn[k] *= n;

        double dqH[4], qMid[4], qNew[4];
        simRotIncrement( vs, 0.5 * dt, dqH );
        simQuatMul( qn,   dqH, qMid );
        simQuatMul( qMid, dqH, qNew );

        // rdot = R(Q) * v_body. R is taken at the MIDDLE of the drift, not at
        // its start: the midpoint rule is self-adjoint and locally O(tau^3),
        // and that is exactly what a symmetric Verlet composition needs from
        // its half-drifts. The pre-step attitude used here before was O(tau^2)
        // and would have held the whole scheme to first order.
        double R[3][3], rd[3];
        quatToR( qMid, R );
        for ( int k = 0 ; k < 3 ; k++ )
            rd[k] = R[k][0]*vs[3] + R[k][1]*vs[4] + R[k][2]*vs[5];

        // Position read before the attitude is written: qOut may alias q.
        double r0 = qs[4], r1 = qs[5], r2 = qs[6];
        for ( int k = 0 ; k < 4 ; k++ ) out[k] = qNew[k];
        out[4] = r0 + rd[0]*dt;
        out[5] = r1 + rd[1]*dt;
        out[6] = r2 + rd[2]*dt;
    }
}

// The two halves of a step, exposed separately because a symmetric (Verlet)
// composition has to evaluate the forces BETWEEN them -- see the header note.
// Both are cheap: the drift is O(nbod) scalar arithmetic with no matrix work
// beyond one 3x3 per floating body, and the kick is an nv-length axpy.
static void simDrift( const Model *mdl, double *q, const double *qd, double tau ) {
    simIntegrateQ( mdl, q, qd, tau, q );
}

static void simKick( const Model *mdl, double *qd, const double *qdd, double dt ) {
    for ( int k = 0 ; k < mdl->nv ; k++ ) qd[k] += qdd[k] * dt;
}

// One semi-implicit Euler step: kick, then drift on the NEW velocity. qdd comes
// from the caller's force model, which is why this takes it rather than calling
// fwd_dyn_crba itself: the contact engine has to run between the kinematics and
// the dynamics. First order, and kept as the control the order gate in
// src/sim/main.c measures the Verlet composition against.
static void simStep( const Model *mdl, double *q, double *qd, const double *qdd, double dt ) {
    simKick( mdl, qd, qdd, dt );
    simDrift( mdl, q, qd, dt );
}

// ---------------------------------------------------------------------------
// Energy, for the conservation checks. Both pieces come straight from the
// spatial quantities denseKinematics already produces.
//
//   KE = 1/2 * sum_i  v_i^T I_i v_i        (v and I both in body i's frame)
//   PE = -sum_i  m_i * g . p_cg_i          (p_cg_i in base coordinates)
//
// This deliberately does NOT go through H: 1/2 qd^T H qd would share the CRBA
// machinery under test, and the point of an energy check is to be independent
// of it.
// ---------------------------------------------------------------------------
// Total linear momentum in base coordinates. Used for the conservation check:
// with contact forces internal and gravity uniform, dP/dt is exactly M*g no
// matter what the contact model does, so any deviation beyond the integrator's
// own O(dt) truncation is a dropped or unbalanced reaction.
static void simMomentum( const Model *mdl, const double *q, const double *qd, double P[3] ) {
    static double Xa[ MODEL_MAX_BODIES ][6][6];
    static double v [ MODEL_MAX_BODIES ][6];
    denseKinematics( mdl, q, qd, Xa, v );

    P[0] = P[1] = P[2] = 0;
    for ( int i = 0 ; i < mdl->nbod ; i++ ) {
        const ModelBody *b = &mdl->body[ i ];
        // Velocity of the CG point, body coordinates: vO + w x cg.
        double vc[3];
        vc[0] = v[i][3] + v[i][1]*b->cg[2] - v[i][2]*b->cg[1];
        vc[1] = v[i][4] + v[i][2]*b->cg[0] - v[i][0]*b->cg[2];
        vc[2] = v[i][5] + v[i][0]*b->cg[1] - v[i][1]*b->cg[0];
        for ( int k = 0 ; k < 3 ; k++ )
            P[k] += b->mass * ( Xa[i][0][k]*vc[0] + Xa[i][1][k]*vc[1] + Xa[i][2][k]*vc[2] );
    }
}

// World position of a POI, from the base->body transforms: p = r + E^T*t.
static void simPoiWorld( const ModelPoi *poi, double Xa[][6][6], double p[3] ) {
    if ( poi->body == -1 ) { for ( int k = 0 ; k < 3 ; k++ ) p[k] = poi->p[k]; return; }
    double r[3];
    xformOrigin( Xa[ poi->body ], r );
    for ( int k = 0 ; k < 3 ; k++ ) {
        double s = r[k];
        for ( int i = 0 ; i < 3 ; i++ ) s += Xa[ poi->body ][i][k] * poi->p[i];
        p[k] = s;
    }
}

// Elastic energy stored in the model's own force elements. Easy to forget and
// impossible to ignore: m113 carries ten torsion bars, six shocks and two idler
// tensioners, and the tensioners alone sit 0.27 m inside their free length. Omit
// this and every joule a spring hands back to the bodies reads as energy
// appearing from nowhere. Dampers store nothing, so only k matters here.
static double simElementEnergy( const Model *mdl, const double *q, double Xa[][6][6] ) {
    double pe = 0;
    for ( int e = 0 ; e < mdl->nelem ; e++ ) {
        const ModelElem *el = &mdl->elem[ e ];
        if ( el->k == 0 ) continue;
        if ( el->type == ELEM_TSDA ) {
            double p1[3], p2[3];
            simPoiWorld( &mdl->poi[ el->poi1 ], Xa, p1 );
            simPoiWorld( &mdl->poi[ el->poi2 ], Xa, p2 );
            double dx = p2[0]-p1[0], dy = p2[1]-p1[1], dz = p2[2]-p1[2];
            double l = sqrt( dx*dx + dy*dy + dz*dz ) - el->ref;
            pe += 0.5 * el->k * l * l;
        } else {
            double dq = q[ mdl->qidx[ el->joint ] ] - el->ref;
            pe += 0.5 * el->k * dq * dq;
        }
    }
    return pe;
}

// Total mechanical energy. keOut, if non-NULL, receives the kinetic part on
// its own -- callers need it as a SCALE, because the absolute value of the
// total is meaningless (the zero of gravitational potential is arbitrary, and
// on a falling vehicle it is also large and drifting).
static double simEnergyK( const Model *mdl, const double *q, const double *qd, double *keOut ) {
    static double Xa[ MODEL_MAX_BODIES ][6][6];   // static: 590 KB, blows the stack
    static double v [ MODEL_MAX_BODIES ][6];
    denseKinematics( mdl, q, qd, Xa, v );

    double ke = 0, pe = simElementEnergy( mdl, q, Xa );
    for ( int i = 0 ; i < mdl->nbod ; i++ ) {
        const ModelBody *b = &mdl->body[ i ];

        double I[6][6], Iv[6];
        mcI( b->mass, b->cg, b->inertia, I );
        m6_mulvec( I, v[i], Iv );
        for ( int k = 0 ; k < 6 ; k++ ) ke += 0.5 * v[i][k] * Iv[k];

        // CG in base coordinates: p = r + E^T * cg, the same recovery
        // poiWorldPos uses (E is Xa's upper-left block, r its frame origin).
        double r[3], p[3];
        xformOrigin( Xa[i], r );
        for ( int k = 0 ; k < 3 ; k++ ) {
            double s = r[k];
            for ( int j = 0 ; j < 3 ; j++ ) s += Xa[i][j][k] * b->cg[j];
            p[k] = s;
        }
        for ( int k = 0 ; k < 3 ; k++ ) pe -= b->mass * mdl->gravity[k] * p[k];
    }
    if ( keOut ) *keOut = ke;
    return ke + pe;
}

static double simEnergy( const Model *mdl, const double *q, const double *qd ) {
    return simEnergyK( mdl, q, qd, NULL );
}
