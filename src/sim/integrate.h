// Plain-C time stepping for the generated forward dynamics.
//
// Until now nothing in build.bat stepped a trajectory: src/validate and
// src/bench evaluate fwd_dyn_crba at random states, and the only integrator in
// the repo lived inside src/chrono/main.cpp, which needs the Chrono SDK to
// compile. The track contact engine cannot be developed or tested without a
// stepper, so integrateQ is ported here with no dependencies beyond the model.
//
// The scheme is semi-implicit (symplectic) Euler, the same one
// stepGenerated() uses in the Chrono cross-check:
//
//     qdd = f( q, qd, ... );   qd += qdd*dt;   q = integrateQ( q, qd, dt );
//
// Velocity is updated first and the NEW velocity advances the position; that
// is what makes the scheme symplectic, so energy oscillates within a bounded
// band instead of growing secularly the way explicit Euler does.
//
// The subtlety this file exists for is integrateQ. For a single-DOF joint q is
// a plain angle or displacement and q += qd*dt is exact. For a FLOATING body
// it is not: qd carries the body's spatial velocity in its OWN frame while q
// carries a base-frame quaternion plus a base-frame position, so the rates are
//
//     Qdot = 0.5 * Q x (0, w_body)        rdot = R(Q) * v_body
//
// Getting this wrong is invisible to src/validate (which never integrates) and
// to the acceleration cross-check (which compares at a fixed state) -- exactly
// the blind spot that once hid the E = R(q) convention bug, so the convention
// note in model.h is the contract here too.

#pragma once
#include <math.h>
#include <string.h>
#include "../model/kinematics.h"

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

        // Quaternion rate Qdot = 0.5 * Q x (0, w_body), with the Hamilton
        // product written out (w,x,y,z order, as everywhere else in the repo).
        double w = qs[0], x = qs[1], y = qs[2], z = qs[3];
        double n = 1.0 / sqrt( w*w + x*x + y*y + z*z );
        w *= n; x *= n; y *= n; z *= n;

        double ox = vs[0], oy = vs[1], oz = vs[2];   // w_body
        double dw = 0.5 * ( -x*ox - y*oy - z*oz );
        double dx = 0.5 * (  w*ox + y*oz - z*oy );
        double dy = 0.5 * (  w*oy + z*ox - x*oz );
        double dz = 0.5 * (  w*oz + x*oy - y*ox );

        double nw = w + dw*dt, nx = x + dx*dt, ny = y + dy*dt, nz = z + dz*dt;
        double nn = 1.0 / sqrt( nw*nw + nx*nx + ny*ny + nz*nz );

        // rdot = R(Q) * v_body, using the PRE-step attitude (the position and
        // attitude updates are both first order, so mixing them buys nothing).
        double R[3][3];
        double qn4[4] = { w, x, y, z };
        quatToR( qn4, R );
        double rd[3];
        for ( int k = 0 ; k < 3 ; k++ )
            rd[k] = R[k][0]*vs[3] + R[k][1]*vs[4] + R[k][2]*vs[5];

        out[0] = nw*nn; out[1] = nx*nn; out[2] = ny*nn; out[3] = nz*nn;
        out[4] = qs[4] + rd[0]*dt;
        out[5] = qs[5] + rd[1]*dt;
        out[6] = qs[6] + rd[2]*dt;
    }
}

// One semi-implicit Euler step. qdd comes from the caller's force model, which
// is why this takes it rather than calling fwd_dyn_crba itself: the contact
// engine has to run between the kinematics and the dynamics.
static void simStep( const Model *mdl, double *q, double *qd, const double *qdd, double dt ) {
    for ( int k = 0 ; k < mdl->nv ; k++ ) qd[k] += qdd[k] * dt;
    simIntegrateQ( mdl, q, qd, dt, q );
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
