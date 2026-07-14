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

    System sys( mdl.nbod );
    for ( int i = 0 ; i < mdl.nbod ; i++ ) {
        ModelBody *b = &mdl.body[ i ];

        // The symbolic pipeline currently only knows how to emit S, xJoint,
        // and vJ for a z-axis revolute joint. The file format already
        // expresses more (prismatic, other axes, cylindrical, spherical);
        // widen this check as the generator learns each one.
        if ( b->jtype != JOINT_REVOLUTE || b->axis != 'z' ) {
            fprintf( stderr, "codegen: body %d: only 'joint revolute z' is implemented so far\n", i );
            return 1;
        }

        sys.parent[ i ]   = (int16_t)b->parent;
        sys.children[ i ] = i + 1;

        string q = "q"s + to_string( i );
        sys.xTree[ i ]  = Xform( format( "Xtree_{}", i ), b->tree[0], b->tree[1], b->tree[2] );
        sys.xJoint[ i ] = Xform( XFORM_RZ, format( "Xjoint_{}", i ), q );

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

    // sys.grav seeds the RNEA acceleration recursion with -a_grav, so negate
    // the physical gravity from the file. numStr keeps zero components "0".
    for ( int k = 0 ; k < 3 ; k++ ) sys.grav.v[ 3 + k ] = numStr( -mdl.gravity[ k ] );

    const char *outPath = argc > 2 ? argv[ 2 ] : "dynamics.c";
    FILE *f = fopen( outPath, "w" );
    if ( !f ) {
        fprintf( stderr, "codegen: cannot write '%s'\n", outPath );
        return 1;
    }
    gen_prelude( f );
    gen_fwd_dyn_crba( f, &sys );
    gen_body_inertias( f, &sys );
    gen_rnea_C( f, &sys );
    gen_crba_H( f, &sys );
    gen_solve_qdd( f, &sys );
    fclose( f );

    return 0;
}
