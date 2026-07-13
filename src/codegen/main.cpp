#include <inttypes.h>
#include <format>
#include "system.h"
#include "vec3.h"
#include "inertia.h"
#include "gen_prelude.h"
#include "gen_fwd_dyn_crba.h"
#include "gen_rnea.h"
#include "gen_crba.h"
#include "gen_solve_qdd.h"

#define NBOD 3

Vec3 cgs[] = {
    { 0.5, 0.0, 0.0 },
    { 1.0, 0.5, 0.0 },
    { 1.5, 1.0, 0.0 },
};

Vec3 jointAxes[] = {
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 0.0, 1.0 },
};

Vec3 jointLocs[] = {
    { 0.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 1.0, 0.0 },
};

double masses[] = {
    1.0,
    1.0,
    1.0
};

Inertia inertias[] = {
    { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 },
    { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 },
    { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 },
};

// Frame frames[ NBOD ] = { 0 };

int main() {
    System sys( NBOD );
    for ( int i = 0 ; i < NBOD ; i++ ) {
        sys.parent[ i ]   = i - 1;
        sys.children[ i ] = i + 1;

        string q = "q"s + to_string( i );

        auto tName = std::format( "Xtree_{}", i );
        auto jName = std::format( "Xjoint_{}", i );
        sys.xTree[ i ]       = Xform( XFORM_TX, tName, "1"s );
        sys.xJoint[ i ]      = Xform( XFORM_RZ, jName, q );

        sys.bodies[ i ].m = masses[ i ];
        sys.bodies[ i ].cg[ 0 ] = cgs[ i ].x;
        sys.bodies[ i ].cg[ 1 ] = cgs[ i ].y;
        sys.bodies[ i ].cg[ 2 ] = cgs[ i ].z;
        sys.bodies[ i ].I = inertias[ i ];
        sys.hasFext[ i ] = 1; // accept an external force on every body; set 0 to omit the code for a body
        // Vec3 newx = v3_cross( &pframe.z, ja );

        // rlink = v3_sub( jp, &plink );
        // plink = *jp;
        // if ( v3_isZero( &rlink ) ) {
        //     sys.xTree[ i ]->type = XFORM_IDENT;
        // }

    }

    // Physical gravity is (0,0,-9.81); sys.grav must hold the negated value
    // used to seed the RNEA acceleration recursion (see System::grav).
    sys.grav.v[ 5 ] = "9.81";

    FILE *f = fopen( "dynamics.c", "w" );
    gen_prelude( f );
    gen_fwd_dyn_crba( f, &sys );
    gen_body_inertias( f, &sys );
    gen_rnea_C( f, &sys );
    gen_crba_H( f, &sys );
    gen_solve_qdd( f, &sys );
    fclose( f );

    return 0;
}