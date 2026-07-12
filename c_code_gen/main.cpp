#include <inttypes.h>
#include <format>
#include "system.h"
#include "vec3.h"
#include "inertia.h"
#include "gen_fwd_dyn_crba.h"

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
        Vec3 rlink;
        Vec3 *jp = &jointLocs[ i ];
        Vec3 *ja = &jointAxes[ i ];
        Vec3 *jpn = &jointLocs[ i+1 ];
        Vec3 *jan = &jointAxes[ i+1 ];

        sys.parent[ i ]   = i - 1;
        sys.children[ i ] = i + 1;

        string q = "q"s + to_string( i );

        auto tName = std::format( "Xtree_{}", i );
        auto jName = std::format( "Xjoint_{}", i );
        sys.xTree[ i ]       = Xform( XFORM_TX, tName, "1"s );
        sys.xJoint[ i ]      = Xform( XFORM_RZ, jName, q );

        // sys.bodies[ i ].m = masses[ i ];
        // sys.bodies[ i ].m = masses[ i ];
        // sys.bodies[ i ].I = inertias[ i ];
        // Vec3 newx = v3_cross( &pframe.z, ja );

        // rlink = v3_sub( jp, &plink );
        // plink = *jp;
        // if ( v3_isZero( &rlink ) ) {
        //     sys.xTree[ i ]->type = XFORM_IDENT;
        // }

    }
    FILE *f = fopen( "dynamics.c", "w" );
    gen_fwd_dyn_crba( f, &sys );
    fclose( f );

    return 0;
}