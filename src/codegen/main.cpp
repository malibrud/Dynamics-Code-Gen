#include <inttypes.h>
#include <format>
#include <stdio.h>
#include "../model/model.h"
#include "system.h"
#include "vec3.h"
#include "inertia.h"
#include "gen_prelude.h"
#include "gen_fwd_dyn_crba.h"
#include "gen_rnea.h"
#include "gen_crba.h"
#include "gen_loops.h"
#include "gen_solve_qdd.h"

int main( int argc, char **argv ) {
    if ( argc < 2 ) {
        fprintf( stderr, "usage: codegen <mechanism-file> [output.c]\n" );
        return 1;
    }

    Model mdl;
    char err[ 256 ];
    if ( modelLoad( argv[ 1 ], &mdl, err, sizeof err ) ) {
        fprintf( stderr, "codegen: %s\n", err );
        return 1;
    }
    const char *unsup = modelUnsupported( &mdl );
    if ( unsup ) {
        fprintf( stderr, "codegen: %s\n", unsup );
        return 1;
    }

    System sys( mdl.nbod );
    for ( int i = 0 ; i < mdl.nbod ; i++ ) {
        ModelBody *b = &mdl.body[ i ];

        sys.parent[ i ]   = (int16_t)b->parent;
        sys.children[ i ] = i + 1;

        // Map the model's joint onto the symbolic joint transform and record
        // type/axis so gen_fwd_dyn_crba can build S and vJ. For a floating
        // joint the transform depends on emitted quaternion temps, so
        // gen_fwd_dyn_crba builds xJoint itself.
        sys.jaxis[ i ] = b->axis;
        sys.xTree[ i ] = Xform( format( "Xtree_{}", i ), b->tree[0], b->tree[1], b->tree[2] );
        if ( b->jtype == JOINT_FLOATING ) {
            sys.jtype[ i ] = J_FLOATING;
        } else {
            Xtype xt;
            if ( b->jtype == JOINT_REVOLUTE ) {
                xt = b->axis == 'x' ? XFORM_RX : b->axis == 'y' ? XFORM_RY : XFORM_RZ;
                sys.jtype[ i ] = J_REVOLUTE;
            } else { // JOINT_PRISMATIC (spherical/cylindrical rejected above)
                xt = b->axis == 'x' ? XFORM_TX : b->axis == 'y' ? XFORM_TY : XFORM_TZ;
                sys.jtype[ i ] = J_PRISMATIC;
            }
            string q = "q"s + to_string( i );
            sys.xJoint[ i ] = Xform( xt, format( "Xjoint_{}", i ), q );
        }
        sys.qidx[ i ] = mdl.qidx[ i ];
        sys.vidx[ i ] = mdl.vidx[ i ];

        sys.bodies[ i ].m = b->mass;
        for ( int k = 0 ; k < 3 ; k++ ) sys.bodies[ i ].cg[ k ] = b->cg[ k ];
        sys.bodies[ i ].I.xx = b->inertia[ 0 ];
        sys.bodies[ i ].I.yy = b->inertia[ 1 ];
        sys.bodies[ i ].I.zz = b->inertia[ 2 ];
        sys.bodies[ i ].I.xy = b->inertia[ 3 ];
        sys.bodies[ i ].I.yz = b->inertia[ 4 ];
        sys.bodies[ i ].I.zx = b->inertia[ 5 ];
        sys.hasFext[ i ] = (uint8_t)b->hasFext;
    }
    sys.nq = mdl.nq;
    sys.nv = mdl.nv;

    // Loop joints (RBDA ch. 8): constraints only, no coordinates.
    sys.stab = mdl.stab;
    for ( int l = 0 ; l < mdl.nloop ; l++ ) {
        ModelLoop *lp = &mdl.loop[ l ];
        System::SysLoop sl;
        sl.pred = lp->pred;
        sl.succ = lp->succ;
        sl.axis = lp->axis;
        switch ( lp->jtype ) {
            case JOINT_REVOLUTE:    sl.jtype = J_REVOLUTE;    break;
            case JOINT_PRISMATIC:   sl.jtype = J_PRISMATIC;   break;
            case JOINT_CYLINDRICAL: sl.jtype = J_CYLINDRICAL; break;
            case JOINT_SPHERICAL:   sl.jtype = J_SPHERICAL;   break;
            case JOINT_DISTANCE:    sl.jtype = J_DISTANCE;    break;
            default:                sl.jtype = J_POINTONLINE; break;
        }
        for ( int k = 0 ; k < 3 ; k++ ) {
            sl.ptree[ k ] = lp->ptree[ k ];
            sl.stree[ k ] = lp->stree[ k ];
        }
        sl.length = lp->length;
        sys.loops.push_back( sl );
    }

    // sys.grav seeds the RNEA acceleration recursion with -a_grav, so negate
    // the physical gravity from the file. numStr keeps zero components "0".
    for ( int k = 0 ; k < 3 ; k++ ) sys.grav.v[ 3 + k ] = numStr( -mdl.gravity[ k ] );

    const char *outPath = argc > 2 ? argv[ 2 ] : "dynamics.c";
    FILE *f = fopen( outPath, "w" );
    if ( !f ) {
        fprintf( stderr, "codegen: cannot write '%s'\n", outPath );
        return 1;
    }
    gen_prelude( f, mdl.nloop > 0 );
    gen_fwd_dyn_crba( f, &sys );
    gen_body_inertias( f, &sys );
    gen_rnea_C( f, &sys );
    gen_crba_H( f, &sys );
    gen_loops( f, &sys );
    gen_solve_qdd( f, &sys );
    fclose( f );

    return 0;
}
