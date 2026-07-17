#pragma once

#include <stdio.h>
#include "system.h"
#include "prln.h"
#include "svec.h"
#include "inertia.h"

// Emits the loop-closure constraint quantities of RBDA ch. 8: the constraint
// matrix K, and the full acceleration-constraint RHS k + kstab (eq 8.35),
// following the algorithm of Table 8.4 in base coordinates. Everything is
// symbolic, so the epsilon structure of K (RBDA eq 8.18-8.20: a loop's row
// block only touches tree joints on the two paths from the loop's root to
// the predecessor and successor) comes out as structural zeros, and the
// planar/axis sparsity of the mechanism folds through exactly like it does
// in the tree dynamics.
//
// Table 8.1's position-error formulae assume the loop joint's axis is z, so
// the declared x/y axes are folded into XP/XS as fixed 90-degree rotations
// (signed permutation matrices -- free after symbolic simplification), and a
// single canonical jcalc serves all axes. This mirrors axisAlignE() in
// src/validate/dense.h.
//
// T is constant in successor coordinates for every Table 8.1 joint type and
// is transformed with Xs^T; point-on-line's constraint directions follow
// the predecessor's axis instead, so it uses Xp^T (its T is constant in
// predecessor coordinates; RBDA eq 8.47 is unchanged either way -- see
// loopJcalc in src/validate/dense.h for the full argument).
//
// distance (a rigid rod between two anchor points, nc = 1) is fixed to
// NEITHER body: its constraint direction rho = p_s - p_p turns with the
// mechanism. It still fits the single-T machinery because the constraint
// force is a pure force along the line through BOTH anchors, so
// T = (p_p x rho; rho) = (p_s x rho; rho) serves predecessor and successor
// alike -- but its acceleration RHS carries the rod's own direction rate
// instead of Table 8.4's vs x vp term (see the branch below).
//
// Must run after gen_rnea_C(): it uses sys->v[] / sys->a[] (the RNEA's
// aJ = 0 acceleration IS Table 8.4's avp; the -a_grav seed offset cancels
// in the successor-minus-predecessor difference, RBDA sec 8.6 footnote) and
// the emitted xFromBase[] transforms.

static int loopNc( uint8_t jtype ) {
    switch ( jtype ) {
        case J_REVOLUTE:
        case J_PRISMATIC:   return 5;
        case J_CYLINDRICAL: return 4;
        case J_SPHERICAL:   return 3;
        case J_POINTONLINE: return 2;
        case J_DISTANCE:    return 1;
        default:            return 0;
    }
}

// Fixed transform from a body frame to a loop joint frame: translation t in
// body coords, plus the axis-alignment rotation putting the joint axis on
// canonical z (E entries are 0/1 permutations, so downstream symbolic math
// stays sparse). Identical to axisAlignE() + plux() in the dense reference.
static Xform loopFrameXform( string name, char axis, const double t[3] ) {
    Xform X( name, t[0], t[1], t[2] ); // identity E, numeric r
    if ( axis == 'x' ) {        // E = R_axis^T, R_axis: jx->by, jy->bz, jz->bx
        X.E[0][0] = "0"; X.E[0][1] = "1"; X.E[0][2] = "0";
        X.E[1][0] = "0"; X.E[1][1] = "0"; X.E[1][2] = "1";
        X.E[2][0] = "1"; X.E[2][1] = "0"; X.E[2][2] = "0";
    } else if ( axis == 'y' ) { // R_axis: jx->bz, jy->bx, jz->by
        X.E[0][0] = "0"; X.E[0][1] = "0"; X.E[0][2] = "1";
        X.E[1][0] = "1"; X.E[1][1] = "0"; X.E[1][2] = "0";
        X.E[2][0] = "0"; X.E[2][1] = "1"; X.E[2][2] = "0";
    }
    return X;
}

// Position error delta = T^T d for the canonical (z-axis) loop joint types,
// RBDA Table 8.1 (X is the measured sXp as a 6x6 of expressions, 0-based
// indices). point-on-line is not in the table; its delta is the exact
// in-plane position of the successor origin, recovered from rho~ = -E^T B
// (see the dense reference). distance never comes here (no loop-joint frame
// at all). Returns nc and fills delta[0..nc-1].
static int loopDelta( uint8_t jtype, array<array<string,6>,6> &X, string delta[6] ) {
    switch ( jtype ) {
    case J_REVOLUTE:
        delta[0] = X[1][2];
        delta[1] = negate( X[0][2] );
        delta[2] = X[4][2];
        delta[3] = negate( X[3][2] );
        delta[4] = addTerms( multFactors( X[3][0], X[1][0] ), multFactors( X[3][1], X[1][1] ) );
        return 5;
    case J_PRISMATIC:
        delta[0] = multFactors( "0.5", subTerms( X[1][2], X[2][1] ) );
        delta[1] = multFactors( "0.5", subTerms( X[2][0], X[0][2] ) );
        delta[2] = multFactors( "0.5", subTerms( X[0][1], X[1][0] ) );
        delta[3] = X[4][2];
        delta[4] = negate( X[3][2] );
        return 5;
    case J_CYLINDRICAL:
        delta[0] = X[1][2];
        delta[1] = negate( X[0][2] );
        delta[2] = X[4][2];
        delta[3] = negate( X[3][2] );
        return 4;
    case J_SPHERICAL:
        delta[0] = addTerms( addTerms( multFactors( X[4][0], X[2][0] ), multFactors( X[4][1], X[2][1] ) ),
                             multFactors( X[4][2], X[2][2] ) );
        delta[1] = addTerms( addTerms( multFactors( X[5][0], X[0][0] ), multFactors( X[5][1], X[0][1] ) ),
                             multFactors( X[5][2], X[0][2] ) );
        delta[2] = addTerms( addTerms( multFactors( X[3][0], X[1][0] ), multFactors( X[3][1], X[1][1] ) ),
                             multFactors( X[3][2], X[1][2] ) );
        return 3;
    default: { // J_POINTONLINE
        // rho_x = -sum_k E[k][2]*B[k][1], rho_y = -sum_k E[k][0]*B[k][2]
        string rx = "0", ry = "0";
        for ( int k = 0 ; k < 3 ; k++ ) {
            rx = subTerms( rx, multFactors( X[k][2], X[k+3][1] ) );
            ry = subTerms( ry, multFactors( X[k][0], X[k+3][2] ) );
        }
        delta[0] = rx;
        delta[1] = ry;
        return 2;
    }
    }
}

// Canonical constraint-force subspace: which rows of the (transposed) frame
// transform become the base-coordinate T columns. A canonical T column e_r
// pulled through X^T is just row r of the 6x6 X, so no multiplication is
// ever emitted. Returns nc and fills rows[0..nc-1].
static int loopTrows( uint8_t jtype, int rows[6] ) {
    switch ( jtype ) {
    case J_REVOLUTE:    rows[0]=0; rows[1]=1; rows[2]=3; rows[3]=4; rows[4]=5; return 5;
    case J_PRISMATIC:   rows[0]=0; rows[1]=1; rows[2]=2; rows[3]=3; rows[4]=4; return 5;
    case J_CYLINDRICAL: rows[0]=0; rows[1]=1; rows[2]=3; rows[3]=4;            return 4;
    case J_SPHERICAL:   rows[0]=3; rows[1]=4; rows[2]=5;                       return 3;
    default:            rows[0]=3; rows[1]=4;                                  return 2; // J_POINTONLINE
    }
}

// Emits each element of a 3-vector of expressions as a named temp, then
// replaces the element with the temp's name (the Svec::emit convention:
// literals 0 and 1 stay as-is to preserve sparsity).
static void emitVec3( string name, string v[3] ) {
    for ( int k = 0 ; k < 3 ; k++ ) {
        prln( "double %s%d = %s;", name.c_str(), k, v[k].c_str() );
        if ( v[k] != "0" && v[k] != "1" ) v[k] = format( "{}{}", name, k );
    }
}

static string dot3( const string a[3], const string b[3] ) {
    string s = "0";
    for ( int k = 0 ; k < 3 ; k++ ) s = addTerms( s, multFactors( a[k], b[k] ) );
    return s;
}

// A loop end's anchor point in base coordinates: p = r + E^T p_local for a
// body end (X = xFromBase maps base->body, so E rotates base coords into
// body coords and r is the body origin in base coords), or the numeric
// anchor itself for a base end.
static void loopAnchorBase( System *sys, int body, const double t[3], string p[3] ) {
    for ( int k = 0 ; k < 3 ; k++ ) {
        if ( body == -1 ) { p[k] = numStr( t[k] ); continue; }
        Xform &X = sys->xFromBase[ body ];
        string s = X.r[k];
        for ( int i = 0 ; i < 3 ; i++ )
            s = addTerms( s, multFactors( X.E[i][k], numStr( t[i] ) ) );
        p[k] = s;
    }
}

void gen_loops( FILE *f, System *sys ) {
    if ( sys->loops.empty() ) return;
    prln_file = f;
    int NL = (int)sys->loops.size();
    int NV = sys->nv;

    sys->nc = 0;
    for ( auto &lp : sys->loops ) sys->nc += loopNc( lp.jtype );
    sys->Kentry.assign( sys->nc, vector<string>( NV, "0"s ) );
    sys->kfulName.assign( sys->nc, ""s );

    string twoAlpha = numStr( 2.0 / sys->stab );
    string betaSq   = numStr( 1.0 / ( sys->stab * sys->stab ) );

    prln( "" );
    prln( "// Loop-closure constraints (RBDA ch. 8, Table 8.4, base coordinates):" );
    prln( "// K rows and the RHS k + kstab of K*qdd = k + kstab (eq 8.35), with" );
    prln( "// Baumgarte stabilization alpha = beta = 1/%g (sec 8.3).", sys->stab );
    prln( "" );

    int row0 = 0;
    for ( int l = 0 ; l < NL ; l++ ) {
        System::SysLoop &lp = sys->loops[ l ];
        int ncl = loopNc( lp.jtype );
        prln( "// Loop %d: predecessor %d, successor %d, %d constraints", l, lp.pred, lp.succ, ncl );

        // Body velocity and velocity-product acceleration in base coords.
        // The fixed base has v = 0 and avp = -a_grav (the RNEA seed): NOT
        // zero, or the -ag offset fails to cancel in as - ap.
        Svec vp( format( "vpL{}_", l ) ), ap = sys->grav, vs( format( "vsL{}_", l ) ), as = sys->grav;
        if ( lp.pred != -1 ) {
            Xform inv = invertXform( format( "XbinvP{}", l ), sys->xFromBase[ lp.pred ] );
            inv.emit();
            vp = fwdXformMotion( format( "vpL{}_", l ), inv, sys->v[ lp.pred ] );
            ap = fwdXformMotion( format( "apL{}_", l ), inv, sys->a[ lp.pred ] );
            vp.emit();
            ap.emit();
        }
        if ( lp.succ != -1 ) {
            Xform inv = invertXform( format( "XbinvS{}", l ), sys->xFromBase[ lp.succ ] );
            inv.emit();
            vs = fwdXformMotion( format( "vsL{}_", l ), inv, sys->v[ lp.succ ] );
            as = fwdXformMotion( format( "asL{}_", l ), inv, sys->a[ lp.succ ] );
            vs.emit();
            as.emit();
        }
        Svec vrel = svecSub( format( "vrel{}_", l ), vs, vp );
        vrel.emit();

        vector<Svec> T( ncl ); // constraint-force columns, base coordinates

        if ( lp.jtype == J_DISTANCE ) {
            // Rod constraint |p_s - p_p| = d, carried as the squared-form
            // position error delta = (rho.rho - d^2)/2 with rho = p_s - p_p:
            // polynomial in the state (no square roots), and its gradient is
            // exactly the K row that T produces, so closure projection works
            // unchanged. T = (p_p x rho; rho) = (p_s x rho; rho) is one
            // spatial force along the rod's line, valid for both bodies.
            // T is fixed to neither body, so instead of Table 8.4's vs x vp
            // term the acceleration RHS is derived directly:
            //   deltaddot = T.(as - ap) + rho.(ws x psdot - wp x ppdot)
            //             + rhodot.rhodot                     (+ K qdd terms)
            // with pdot = v0 + w x p the anchor point velocities.
            string pp[3], ps[3], rho[3];
            loopAnchorBase( sys, lp.pred, lp.ptree, pp );
            loopAnchorBase( sys, lp.succ, lp.stree, ps );
            emitVec3( format( "ppL{}_", l ), pp );
            emitVec3( format( "psL{}_", l ), ps );
            for ( int k = 0 ; k < 3 ; k++ ) rho[k] = subTerms( ps[k], pp[k] );
            emitVec3( format( "rhoL{}_", l ), rho );

            // Prefer the base-side anchor for the moment arm (numeric
            // components); either anchor gives the same T.
            const string *anchor = ( lp.pred != -1 && lp.succ == -1 ) ? ps : pp;
            string ang[3];
            cross3( anchor, rho, ang );
            T[0] = Svec( format( "TL{}_0_", l ) );
            for ( int k = 0 ; k < 3 ; k++ ) { T[0].v[k] = ang[k]; T[0].v[3+k] = rho[k]; }
            T[0].emit();

            // Anchor point velocities and the rod direction rate.
            string ppd[3], psd[3], rhod[3], wxp[3];
            cross3( vp.v, pp, wxp );
            for ( int k = 0 ; k < 3 ; k++ ) ppd[k] = addTerms( vp.v[3+k], wxp[k] );
            emitVec3( format( "ppdL{}_", l ), ppd );
            cross3( vs.v, ps, wxp );
            for ( int k = 0 ; k < 3 ; k++ ) psd[k] = addTerms( vs.v[3+k], wxp[k] );
            emitVec3( format( "psdL{}_", l ), psd );
            for ( int k = 0 ; k < 3 ; k++ ) rhod[k] = subTerms( psd[k], ppd[k] );
            emitVec3( format( "rhodL{}_", l ), rhod );

            string wxpsd[3], wxppd[3];
            cross3( vs.v, psd, wxpsd );
            cross3( vp.v, ppd, wxppd );
            string kvel = dot3( rhod, rhod );
            for ( int k = 0 ; k < 3 ; k++ )
                kvel = addTerms( kvel, multFactors( rho[k], subTerms( wxpsd[k], wxppd[k] ) ) );

            Svec dAvp = svecSub( format( "dAvp{}_", l ), as, ap );
            dAvp.emit();

            string delta = multFactors( "0.5", subTerms( dot3( rho, rho ),
                                                         numStr( lp.length * lp.length ) ) );
            string kstab = addTerms( multFactors( twoAlpha, svecDot( T[0], vrel ) ),
                                     multFactors( betaSq, delta ) );
            string kful = subTerms( negate( addTerms( svecDot( T[0], dAvp ), kvel ) ), kstab );
            string name = format( "kful_{}", row0 );
            prln( "double %s = %s;", name.c_str(), kful.c_str() );
            sys->kfulName[ row0 ] = name;
        } else {
            // Xp = XP * pred_X_base, Xs = XS * succ_X_base (Table 8.4's Xp/Xs).
            Xform XP = loopFrameXform( format( "XPf{}", l ), lp.axis, lp.ptree );
            Xform XS = loopFrameXform( format( "XSf{}", l ), lp.axis, lp.stree );
            Xform Xp = lp.pred == -1 ? XP : compose( format( "Xp{}", l ), XP, sys->xFromBase[ lp.pred ] );
            Xform Xs = lp.succ == -1 ? XS : compose( format( "Xs{}", l ), XS, sys->xFromBase[ lp.succ ] );
            Xp.baseName = format( "Xp{}", l );
            Xs.baseName = format( "Xs{}", l );
            Xp.emit();
            Xs.emit();

            // sXp = Xs * Xp^-1: the measured predecessor->successor transform of
            // the loop joint, input to the Table 8.1 position-error formulae.
            Xform XpInv = invertXform( format( "XpInv{}", l ), Xp );
            XpInv.emit();
            Xform sXp = compose( format( "sXp{}", l ), Xs, XpInv );
            sXp.emit();
            auto sXp6 = xformTo6x6( sXp );

            string delta[6];
            loopDelta( lp.jtype, sXp6, delta );

            // T columns in base coordinates: rows of the 6x6 of Xs (successor-
            // fixed types) or Xp (point-on-line; predecessor-fixed).
            int rows[6];
            loopTrows( lp.jtype, rows );
            Xform &Tsrc = lp.jtype == J_POINTONLINE ? Xp : Xs;
            auto Tsrc6 = xformTo6x6( Tsrc );
            for ( int c = 0 ; c < ncl ; c++ ) {
                T[c] = Svec( format( "TL{}_{}_", l, c ) );
                for ( int r = 0 ; r < 6 ; r++ ) T[c].v[r] = Tsrc6[ rows[c] ][ r ];
                T[c].emit();
            }

            // kful = -T^T (as - ap + vs x vp) - ( 2*alpha*T^T(vs - vp) + beta^2*delta ),
            // which is k + kstab of eq 8.35 (Table 8.4's local kstab is the
            // negative of eq 8.22's, and the table subtracts it).
            Svec asMinusAp = svecSub( format( "dAvp{}_", l ), as, ap );
            Svec vsxvp = svCrossMotion( format( "vsxvp{}_", l ), vs, vp );
            Svec resid = svecAdd( format( "resid{}_", l ), asMinusAp, vsxvp );
            resid.emit();

            for ( int c = 0 ; c < ncl ; c++ ) {
                string kstab = addTerms( multFactors( twoAlpha, svecDot( T[c], vrel ) ),
                                         multFactors( betaSq, delta[c] ) );
                string kful = subTerms( negate( svecDot( T[c], resid ) ), kstab );
                string name = format( "kful_{}", row0 + c );
                prln( "double %s = %s;", name.c_str(), kful.c_str() );
                sys->kfulName[ row0 + c ] = name;
            }
        }

        // K row block: walk both paths to the common ancestor (parents have
        // smaller indices; the base is -1, so stepping the larger index meets
        // at the nearest common ancestor). -T^T 0X_i S_i on the predecessor
        // side, + on the successor side.
        int i = lp.pred, j = lp.succ;
        while ( i != j ) {
            bool onPred = i > j;
            int b = onPred ? i : j;

            Xform inv = invertXform( format( "XbinvK{}_{}", l, b ), sys->xFromBase[ b ] );
            inv.emit();

            // Base-coordinate motion-subspace columns of tree joint b.
            int ndof = sys->jtype[ b ] == J_FLOATING ? 6 : 1;
            for ( int col = 0 ; col < ndof ; col++ ) {
                Svec Scol( format( "S0L{}_{}_{}_", l, b, col ) );
                if ( sys->jtype[ b ] == J_FLOATING ) {
                    Svec e( ""s );
                    e.v[ col ] = "1";
                    Scol = fwdXformMotion( Scol.baseName, inv, e );
                } else {
                    Scol = fwdXformMotion( Scol.baseName, inv, sys->S[ b ] );
                }
                Scol.emit();
                for ( int c = 0 ; c < ncl ; c++ ) {
                    string e = svecDot( T[c], Scol );
                    if ( onPred ) e = negate( e );
                    if ( e == "0" ) continue; // structural zero stays out of the array
                    string name = format( "K_{}_{}", row0 + c, sys->vidx[ b ] + col );
                    prln( "double %s = %s;", name.c_str(), e.c_str() );
                    sys->Kentry[ row0 + c ][ sys->vidx[ b ] + col ] = name;
                }
            }
            if ( onPred ) i = sys->parent[ i ];
            else          j = sys->parent[ j ];
        }
        prln( "" );
        row0 += ncl;
    }
}
