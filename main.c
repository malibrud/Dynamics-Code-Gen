#include <inttypes.h>
#include "system.h"
#include "fwd_dyn_crba.h"

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

Frame frames[ NBOD ] = { 0 };

int main() {
    System sys = { 0 };
    sys_init( &sys, 3 );
    x_init( sys.arena );

    Vec3 plink = { 0 };
    Frame pframe = {
        { 0, 0, 0 },
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 },
    };

    for ( int16_t i = 0 ; i < NBOD ; i++ ) {
        Vec3 rlink;
        Vec3 *jp = &jointLocs[ i ];
        Vec3 *ja = &jointAxes[ i ];
        Vec3 *jpn = &jointLocs[ i+1 ];
        Vec3 *jan = &jointAxes[ i+1 ];

        sys.parent[ i ]   = i - 1;
        sys.children[ i ] = i + 1;

        sys.xTree[ i ]       = x_newTx( 1.0 );
        sys.xJoint[ i ]      = x_newRz( 0.0 );
        sys.xFromParent[ i ] = x_newTxRz( 1, 0 );
        sys.xFromBase[ i ]   = x_newPxy( 1, 0, 0 );

        sys.bodies[ i ].m = masses[ i ];
        sys.bodies[ i ].m = masses[ i ];
        sys.bodies[ i ].I = inertias[ i ];
        Vec3 newx = v3_cross( &pframe.z, ja );

        rlink = v3_sub( jp, &plink );
        plink = *jp;
        if ( v3_isZero( &rlink ) ) {
            sys.xTree[ i ]->type = XFORM_IDENT;
        }

    }

    return 0;
}