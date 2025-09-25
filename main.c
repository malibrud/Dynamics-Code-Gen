#include <inttypes.h>
#include "system.h"

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
    System sys = {};
    sys_init( &sys, 3 );

    Vec3 plink = { 0 };
    Frame pframe = {
        { 0, 0, 0 },
        { 1, 0, 0 },
        { 0, 1, 0 },
        { 0, 0, 1 },
    };

    for ( int i = 0 ; i < NBOD ; i++ ) {
        Vec3 rlink;
        Vec3 *jp = &jointLocs[ i ];
        Vec3 *ja = &jointAxes[ i ];
        Vec3 *jpn = &jointLocs[ i+1 ];
        Vec3 *jan = &jointAxes[ i+1 ];

        if ( i == 0 ) {
            sys.xTree[ i ] = sys_newXform( &sys, XFORM_IDENT );
            

        }
        else {
            Xform_AxisTrans *X = sys_newXform( &sys, XFORM_TX );
            X->d = v3_dist( jp, jpn );
            sys.xTree[ i ] = X;
        }
        sys.xJoint[ i ]= sys_newXform( &sys, XFORM_RZ );

        sys.bodies[ i ].m = masses[ i ];
        sys.bodies[ i ].m = masses[ i ];
        sys.bodies[ i ].I = inertias[ i ];
        Vec3 newx = v3_cross( &pframe.z, ja );

        rlink = v3_sub( jp, &plink );
        plink = *jp;
        if ( v3_isZero( &rlink ) ) {
            sys.xTree[ i ].type = XFORM_IDENT;
        }

    }

    return 0;
}