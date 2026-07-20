#pragma once

#include <stdio.h>
#include "system.h"
#include "prln.h"
#include "svec.h"
#include "inertia.h"
#include "gen_loops.h"

// Shared TSDA kinematics -- anchor world positions, rho, l, 1/l, ldot --
// needed by both gen_elems() (which continues into the force law and
// path-walk) and gen_observe() (stage 5, which only needs l/ldot for the
// observation array). Must run after gen_rnea_C() (needs sys->v[] and
// sys->xFromBase[]); emits into whatever function body is currently open, so
// it works identically in either generated function's own scope (temp names
// may repeat safely across functions). Fills pf3/pt3/rho with the emitted
// anchor/rho temp names and lname/linvname/ldotname with the l/1/l/ldot temp
// names.
static void gen_tsda_kinematics( FILE *f, System *sys, int e,
                                  string pf3[3], string pt3[3], string rho[3],
                                  string &lname, string &linvname, string &ldotname ) {
    prln_file = f;
    System::SysElem &el = sys->elems[ e ];
    System::SysPoi &pf = sys->pois[ el.poi1 ]; // "from" plays the distance
    System::SysPoi &pt = sys->pois[ el.poi2 ]; // loop's predecessor role, "to" the successor role

    loopAnchorBase( sys, pf.body, pf.p, pf3 );
    loopAnchorBase( sys, pt.body, pt.p, pt3 );
    emitVec3( format( "elemPf{}_", e ), pf3 );
    emitVec3( format( "elemPt{}_", e ), pt3 );
    for ( int k = 0 ; k < 3 ; k++ ) rho[k] = subTerms( pt3[k], pf3[k] );
    emitVec3( format( "elemRho{}_", e ), rho );

    string l2 = dot3( rho, rho );
    prln( "double elemL2_%d = %s;", e, l2.c_str() );
    lname = format( "elemL_{}", e );
    prln( "double %s = sqrt( elemL2_%d );", lname.c_str(), e );
    linvname = format( "elemLinv_{}", e );
    prln( "double %s = 1.0 / %s;", linvname.c_str(), lname.c_str() );

    // Anchor body velocities in base coordinates (no acceleration term:
    // elements carry no constraint, so no Baumgarte/avp machinery).
    Svec vp( format( "elemVp{}_", e ) ), vs( format( "elemVs{}_", e ) );
    if ( pf.body != -1 ) {
        Xform inv = invertXform( format( "elemXbinvF{}", e ), sys->xFromBase[ pf.body ] );
        inv.emit();
        vp = fwdXformMotion( format( "elemVp{}_", e ), inv, sys->v[ pf.body ] );
        vp.emit();
    }
    if ( pt.body != -1 ) {
        Xform inv = invertXform( format( "elemXbinvT{}", e ), sys->xFromBase[ pt.body ] );
        inv.emit();
        vs = fwdXformMotion( format( "elemVs{}_", e ), inv, sys->v[ pt.body ] );
        vs.emit();
    }

    // Anchor point velocities (pdot = v0 + w x p; the spatial-velocity
    // convention where p is measured from the BASE origin makes this valid
    // for any point, not just the body's own origin) and the rod direction
    // rate.
    string wxp[3], pfd[3], ptd[3], rhod[3];
    cross3( vp.v, pf3, wxp );
    for ( int k = 0 ; k < 3 ; k++ ) pfd[k] = addTerms( vp.v[3+k], wxp[k] );
    emitVec3( format( "elemPfd{}_", e ), pfd );
    cross3( vs.v, pt3, wxp );
    for ( int k = 0 ; k < 3 ; k++ ) ptd[k] = addTerms( vs.v[3+k], wxp[k] );
    emitVec3( format( "elemPtd{}_", e ), ptd );
    for ( int k = 0 ; k < 3 ; k++ ) rhod[k] = subTerms( ptd[k], pfd[k] );
    emitVec3( format( "elemRhod{}_", e ), rhod );

    ldotname = format( "elemLdot_{}", e );
    prln( "double %s = %s;", ldotname.c_str(), multFactors( dot3( rho, rhod ), linvname ).c_str() );
}

// Emits the generalized force Q contributed by each force element (TSDA/RSDA,
// src/model/model.h), accumulating into sys->Qentry indexed by velocity
// coordinate -- gen_solve_qdd() adds "+ Q_r" to the rhs wherever an entry is
// non-"0". Must run after gen_rnea_C() (needs sys->v[] and sys->xFromBase[])
// and before gen_solve_qdd(); order relative to gen_crba_H()/gen_loops() does
// not matter, since elements share no state with either.
//
// RSDA needs no geometry: it is a torque straight into its joint's own
// coordinate, tau_e = -k*(q - q0) - c*qd + act[slot].
//
// TSDA is the distance loop joint's geometry (gen_loops.h's J_DISTANCE
// branch) with a known force instead of an unknown constraint force: T =
// (p_from x rho; rho) is one spatial force serving both anchors, walked to
// the tree root with sign - on the "from" (predecessor) path and + on the
// "to" (successor) path (RBDA eq 8.18-8.20's epsilon structure), producing
// row[r] = l * dl/dq_r. The generalized force of tension f is then
// Q[r] = -(f/l) * row[r], f = k*(l - l0) + c*ldot + act[slot] -- exact, not
// small-angle.
void gen_elems( FILE *f, System *sys ) {
    if ( sys->elems.empty() ) return;
    prln_file = f;
    int NV = sys->nv;

    sys->Qentry.assign( NV, "0"s );

    prln( "" );
    prln( "// Force elements (TSDA/RSDA): generalized force Q, applied via the same" );
    prln( "// path-walk projection as a distance loop joint (TSDA) or directly into" );
    prln( "// the owning joint's coordinate (RSDA). See src/model/model.h for the" );
    prln( "// force-law sign conventions." );
    prln( "" );

    for ( int e = 0 ; e < (int)sys->elems.size() ; e++ ) {
        System::SysElem &el = sys->elems[ e ];

        if ( el.type == System::SysElem::RSDA ) {
            int b = el.joint;
            string qVar  = format( "q{}", b );
            string qdVar = format( "qd{}", b );

            string term = "0";
            if ( el.k != 0 ) term = subTerms( term, multFactors( numStr( el.k ), subTerms( qVar, numStr( el.ref ) ) ) );
            if ( el.c != 0 ) term = subTerms( term, multFactors( numStr( el.c ), qdVar ) );
            if ( el.hasInput ) term = addTerms( term, format( "act[ {} ]", el.actIdx ) );

            int r = sys->vidx[ b ];
            sys->Qentry[ r ] = addTerms( sys->Qentry[ r ], term );
            continue;
        }

        System::SysPoi &pf = sys->pois[ el.poi1 ];
        System::SysPoi &pt = sys->pois[ el.poi2 ];

        string pf3[3], pt3[3], rho[3], lname, linvname, ldotname;
        gen_tsda_kinematics( f, sys, e, pf3, pt3, rho, lname, linvname, ldotname );

        // Tension law: f = k*(l - l0) + c*ldot + act[slot], positive f pulls
        // the anchors together.
        string fterm = "0";
        if ( el.k != 0 ) fterm = addTerms( fterm, multFactors( numStr( el.k ), subTerms( lname, numStr( el.ref ) ) ) );
        if ( el.c != 0 ) fterm = addTerms( fterm, multFactors( numStr( el.c ), ldotname ) );
        if ( el.hasInput ) fterm = addTerms( fterm, format( "act[ {} ]", el.actIdx ) );
        string fname = format( "elemF_{}", e );
        prln( "double %s = %s;", fname.c_str(), fterm.c_str() );

        string negFOverLName = format( "elemNegFOverL_{}", e );
        prln( "double %s = %s;", negFOverLName.c_str(), negate( multFactors( fname, linvname ) ).c_str() );

        // Prefer the base-side anchor for the moment arm (numeric
        // components); either anchor gives the same T.
        const string *anchor = ( pf.body != -1 && pt.body == -1 ) ? pt3 : pf3;
        string ang[3];
        cross3( anchor, rho, ang );
        Svec T( format( "elemT{}_", e ) );
        for ( int k = 0 ; k < 3 ; k++ ) { T.v[k] = ang[k]; T.v[3+k] = rho[k]; }
        T.emit();

        // Path walk to the tree root (RBDA eq 8.18-8.20's epsilon structure,
        // identical to gen_loops.h's K row block): - on the "from"
        // (predecessor) side, + on the "to" (successor) side.
        int i = pf.body, j = pt.body;
        while ( i != j ) {
            bool onPred = i > j;
            int b = onPred ? i : j;

            Xform inv = invertXform( format( "elemXbinvK{}_{}", e, b ), sys->xFromBase[ b ] );
            inv.emit();

            int ndof = sys->jtype[ b ] == J_FLOATING ? 6 : 1;
            for ( int col = 0 ; col < ndof ; col++ ) {
                Svec Scol( format( "elemS0{}_{}_{}_", e, b, col ) );
                if ( sys->jtype[ b ] == J_FLOATING ) {
                    Svec unit( ""s );
                    unit.v[ col ] = "1";
                    Scol = fwdXformMotion( Scol.baseName, inv, unit );
                } else {
                    Scol = fwdXformMotion( Scol.baseName, inv, sys->S[ b ] );
                }
                Scol.emit();

                string row = svecDot( T, Scol );
                if ( onPred ) row = negate( row );
                if ( row == "0" ) continue;

                int r = sys->vidx[ b ] + col;
                sys->Qentry[ r ] = addTerms( sys->Qentry[ r ], multFactors( negFOverLName, row ) );
            }
            if ( onPred ) i = sys->parent[ i ];
            else          j = sys->parent[ j ];
        }
        prln( "" );
    }

    // Materialize the accumulated Q_r temps -- only for coordinates some
    // element's path touches (structural zeros stay out, same discipline as
    // Kentry).
    for ( int r = 0 ; r < NV ; r++ ) {
        if ( sys->Qentry[ r ] == "0" ) continue;
        string name = format( "Q_{}", r );
        prln( "double %s = %s;", name.c_str(), sys->Qentry[ r ].c_str() );
        sys->Qentry[ r ] = name;
    }
    prln( "" );
}
