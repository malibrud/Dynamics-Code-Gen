#include <inttypes.h>
#include "system.h"

#define NBOD 3

double cgs[][ 3 ] = {
    { 0.5, 0.0, 0.0 },
    { 1.0, 0.5, 0.0 },
    { 1.5, 1.0, 0.0 },
};

double joints[][ 3 ] = {
    { 0.0, 0.0, 0.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 1.0, 0.0 },
};

double masses[] = {
    1.0,
    1.0,
    1.0
};

double inertias[][ 6 ] = {
    { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 },
    { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 },
    { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 },
};

int main() {
    System sys = {};
    sys_init( &sys, 3 );

    for ( int i = 0 ; i < NBOD ; i++ ) {

    }

    return 0;
}