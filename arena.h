#pragma once
#include <stdlib.h>
#include "check.h"

#define ALIGNMENT ( sizeof( void* ) )
#define ALIGN(x) ( ( ( nb + ( ALIGNMENT-1 ) ) / ALIGNMENT ) * ALIGNMENT )
#define ALIGNTO(x,a) ( ( ( nb + ( a-1 ) ) / a ) * a )

typedef struct {
    char *memory;
    size_t capacity;
    size_t size;
} Arena;

Arena *ar_create( size_t initCapacity ) {
    Arena *ar = malloc( sizeof( Arena ) );
    check( ar, "Could not allocate memory for arena." );

    ar->memory = malloc( initCapacity );
    check( ar->memory, "Could not allocate arena buffer." );
    ar->capacity = initCapacity;
    ar->size = 0;
}

void *ar_alloc_aligned( Arena *ar, size_t size, size_t align ) {
    assert( ar );
    assert( size  > 0 );
    assert( align > 0 );

    size_t start = ALIGNTO( ar->size, align );
    check( start + size > ar->capacity, "Arena is out of capacity." );
    ar->size = start + size;
    return ar->memory + start;
}

void *ar_alloc_array( Arena *ar, size_t element_count, size_t element_size ) {
    assert( ar );
    assert( element_count > 0 );
    assert( element_size  > 0 );

    size_t array_size = element_count * element_size;
    return( ar_alloc_aligned( ar, array_size, element_size ) );
}