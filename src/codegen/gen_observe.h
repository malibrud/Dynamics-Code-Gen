#pragma once

#include <stdio.h>
#include "system.h"
#include "prln.h"
#include "svec.h"
#include "xform.h"
#include "gen_kinematics.h"
#include "gen_elems.h" // gen_tsda_kinematics, loopAnchorBase/emitVec3/cross3 (via gen_loops.h)

// Emits fwd_dyn_find( kind, name ): plain strcmp lookup over per-kind const
// name tables. "poi"/"tsda"/"rsda" return the element's out_obs offset,
// "act" returns the act[] slot (only elements declared 'input' get a name
// here). Runtime lookup rather than a generated header, because the eventual
// JIT consumer is a fixed host program that cannot recompile per model.
// Returns -1 if kind or name is not found. Always emitted (even for a model
// with nothing declared), so the export is a stable part of the API.
static void gen_find( FILE *f, System *sys ) {
    prln_file = f;
    prln( "" );

    int npoi = (int)sys->pois.size();
    int ntsda = 0, nrsda = 0, nactNamed = 0;
    for ( auto &el : sys->elems ) {
        if ( el.type == System::SysElem::TSDA ) ntsda++; else nrsda++;
        if ( el.hasInput ) nactNamed++;
    }

    if ( npoi > 0 ) {
        prln( "static const char *fwd_dyn_poi_names[] = {" );
        prln_indent_level++;
        for ( auto &p : sys->pois ) prln( "\"%s\",", p.name.c_str() );
        prln_indent_level--;
        prln( "};" );
        prln( "static const int fwd_dyn_poi_offsets[] = {" );
        prln_indent_level++;
        for ( auto &p : sys->pois ) prln( "%d,", p.obsIdx );
        prln_indent_level--;
        prln( "};" );
    }
    if ( ntsda > 0 ) {
        prln( "static const char *fwd_dyn_tsda_names[] = {" );
        prln_indent_level++;
        for ( auto &el : sys->elems ) if ( el.type == System::SysElem::TSDA ) prln( "\"%s\",", el.name.c_str() );
        prln_indent_level--;
        prln( "};" );
        prln( "static const int fwd_dyn_tsda_offsets[] = {" );
        prln_indent_level++;
        for ( auto &el : sys->elems ) if ( el.type == System::SysElem::TSDA ) prln( "%d,", el.obsIdx );
        prln_indent_level--;
        prln( "};" );
    }
    if ( nrsda > 0 ) {
        prln( "static const char *fwd_dyn_rsda_names[] = {" );
        prln_indent_level++;
        for ( auto &el : sys->elems ) if ( el.type == System::SysElem::RSDA ) prln( "\"%s\",", el.name.c_str() );
        prln_indent_level--;
        prln( "};" );
        prln( "static const int fwd_dyn_rsda_offsets[] = {" );
        prln_indent_level++;
        for ( auto &el : sys->elems ) if ( el.type == System::SysElem::RSDA ) prln( "%d,", el.obsIdx );
        prln_indent_level--;
        prln( "};" );
    }
    if ( nactNamed > 0 ) {
        prln( "static const char *fwd_dyn_act_names[] = {" );
        prln_indent_level++;
        for ( auto &el : sys->elems ) if ( el.hasInput ) prln( "\"%s\",", el.name.c_str() );
        prln_indent_level--;
        prln( "};" );
        prln( "static const int fwd_dyn_act_slots[] = {" );
        prln_indent_level++;
        for ( auto &el : sys->elems ) if ( el.hasInput ) prln( "%d,", el.actIdx );
        prln_indent_level--;
        prln( "};" );
    }

    prln( "" );
    prln( "// Runtime name lookup. kind is \"poi\" | \"tsda\" | \"rsda\" (returns the" );
    prln( "// element's out_obs offset) or \"act\" (returns the act[] slot). Returns" );
    prln( "// -1 if kind or name is not found." );
    prln( "int fwd_dyn_find( const char *kind, const char *name ) {" );
    prln_indent_level++;

    prln( "if ( !strcmp( kind, \"poi\" ) ) {" );
    prln_indent_level++;
    if ( npoi > 0 )
        prln( "for ( int i = 0 ; i < %d ; i++ ) if ( !strcmp( fwd_dyn_poi_names[ i ], name ) ) return fwd_dyn_poi_offsets[ i ];", npoi );
    prln( "return -1;" );
    prln_indent_level--;
    prln( "}" );

    prln( "if ( !strcmp( kind, \"tsda\" ) ) {" );
    prln_indent_level++;
    if ( ntsda > 0 )
        prln( "for ( int i = 0 ; i < %d ; i++ ) if ( !strcmp( fwd_dyn_tsda_names[ i ], name ) ) return fwd_dyn_tsda_offsets[ i ];", ntsda );
    prln( "return -1;" );
    prln_indent_level--;
    prln( "}" );

    prln( "if ( !strcmp( kind, \"rsda\" ) ) {" );
    prln_indent_level++;
    if ( nrsda > 0 )
        prln( "for ( int i = 0 ; i < %d ; i++ ) if ( !strcmp( fwd_dyn_rsda_names[ i ], name ) ) return fwd_dyn_rsda_offsets[ i ];", nrsda );
    prln( "return -1;" );
    prln_indent_level--;
    prln( "}" );

    prln( "if ( !strcmp( kind, \"act\" ) ) {" );
    prln_indent_level++;
    if ( nactNamed > 0 )
        prln( "for ( int i = 0 ; i < %d ; i++ ) if ( !strcmp( fwd_dyn_act_names[ i ], name ) ) return fwd_dyn_act_slots[ i ];", nactNamed );
    prln( "return -1;" );
    prln_indent_level--;
    prln( "}" );

    prln( "return -1;" );
    prln_indent_level--;
    prln( "}" );
}

// Emits fwd_dyn_observe( q, qd, out_obs ): world kinematics of every POI
// (position + velocity) and each force element's state (TSDA l/ldot; RSDA
// angle/rate), at the obsIdx offsets src/model/model.h documents (all POIs,
// then all TSDAs, then all RSDAs), plus fwd_dyn_observe_len and
// fwd_dyn_find(). Runs the shared kinematics + velocity emission
// (gen_kinematics.h, stage 2) in a SEPARATE function scope from
// fwd_dyn_crba -- temp names may repeat safely -- so it builds the same
// state without duplicating any dynamics (CRBA H, RNEA bias force). Skips
// that recompute entirely for a model with nothing declared, so
// element-free models pay nothing for this export beyond the two constants.
// Must run after fwd_dyn_crba's function has been closed (gen_solve_qdd).
void gen_observe( FILE *f, System *sys ) {
    prln_file = f;
    int NB = sys->nBod;
    int obsLen = 6 * (int)sys->pois.size() + 2 * (int)sys->elems.size();

    prln( "" );
    prln( "// World kinematics of every POI (position + velocity) and each force" );
    prln( "// element's state (TSDA l/ldot; RSDA angle/rate), laid out at the obsIdx" );
    prln( "// offsets documented in src/model/model.h: all POIs (6 doubles each: px" );
    prln( "// py pz vx vy vz), then all TSDAs (2 each: l ldot), then all RSDAs (2" );
    prln( "// each: angle rate)." );
    prln( "const int fwd_dyn_observe_len = %d;", obsLen );
    prln( "" );
    prln( "void fwd_dyn_observe( double *q, double *qd, double *out_obs ) {" );
    prln_indent_level++;

    if ( !sys->pois.empty() || !sys->elems.empty() ) {
        gen_kinematics( f, sys );
        for ( int i = 0 ; i < NB ; i++ ) gen_velocity_body( f, sys, i );

        for ( int p = 0 ; p < (int)sys->pois.size() ; p++ ) {
            System::SysPoi &poi = sys->pois[ p ];

            string pos[3];
            loopAnchorBase( sys, poi.body, poi.p, pos );
            emitVec3( format( "poiP{}_", p ), pos );

            Svec v0( format( "poiV{}_", p ) );
            if ( poi.body != -1 ) {
                Xform inv = invertXform( format( "poiXbinv{}", p ), sys->xFromBase[ poi.body ] );
                inv.emit();
                v0 = fwdXformMotion( format( "poiV{}_", p ), inv, sys->v[ poi.body ] );
                v0.emit();
            }
            // World velocity: v0 + w x p (spatial-velocity convention, valid
            // for any point measured from the base origin).
            string wxp[3], vel[3];
            cross3( v0.v, pos, wxp );
            for ( int k = 0 ; k < 3 ; k++ ) vel[k] = addTerms( v0.v[3+k], wxp[k] );
            emitVec3( format( "poiVel{}_", p ), vel );

            for ( int k = 0 ; k < 3 ; k++ ) prln( "out_obs[ %d ] = %s;", poi.obsIdx + k, pos[k].c_str() );
            for ( int k = 0 ; k < 3 ; k++ ) prln( "out_obs[ %d ] = %s;", poi.obsIdx + 3 + k, vel[k].c_str() );
        }

        for ( int e = 0 ; e < (int)sys->elems.size() ; e++ ) {
            System::SysElem &el = sys->elems[ e ];
            if ( el.type == System::SysElem::RSDA ) {
                prln( "out_obs[ %d ] = q%d;", el.obsIdx, el.joint );
                prln( "out_obs[ %d ] = qd%d;", el.obsIdx + 1, el.joint );
                continue;
            }
            string pf3[3], pt3[3], rho[3], lname, linvname, ldotname;
            gen_tsda_kinematics( f, sys, e, pf3, pt3, rho, lname, linvname, ldotname );
            prln( "out_obs[ %d ] = %s;", el.obsIdx, lname.c_str() );
            prln( "out_obs[ %d ] = %s;", el.obsIdx + 1, ldotname.c_str() );
        }
    }

    prln_indent_level--;
    prln( "}" );

    gen_find( f, sys );
}
