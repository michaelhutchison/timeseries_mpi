#pragma once

#include "vec3.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <stdlib.h>
#include <time.h>
#include <cmath>

#include "Object.h"

#define RANDOMDOUBLE(min,max) (  ((double)rand() / RAND_MAX) * (max-min) + min  ) 
#ifndef PI
#define PI 3.14159265
#endif

class Slice_serial;

struct Bounds {
    double nx; // negative x boundary
    double px; // positive x boundary
    double ny; // negative y boundary
    double py; // positive y boundary
    double nz; // negative z boundary
    double pz; // positive z boundary
};
struct Neighbors_rank {
    // Holds the rank of the neighbor at each boundary.
    // Value of -1 indicates no neighbor on that boundary.
    int nx; // negative x neighbor
    int px; // positive x neighbor
    int ny; // negative y neighbor
    int py; // positive y neighbor
    int nz; // negative z neighbor
    int pz; // positive z neighbor
};
struct Neighbors_pointer {
    // Holds the rank of the neighbor at each boundary.
    // Value of 0 indicates no neighbor on that boundary.
    Slice_serial * nx; // negative x neighbor
    Slice_serial * px; // positive x neighbor
    Slice_serial * ny; // negative y neighbor
    Slice_serial * py; // positive y neighbor
    Slice_serial * nz; // negative z neighbor
    Slice_serial * pz; // positive z neighbor
};


class Slice_serial {
public:
    Slice_serial(unsigned r, unsigned s, std::ofstream * fh, double worldSize[3], unsigned sliceRows[3]);
    ~Slice_serial();
    void setTotalObjects(unsigned n) {nTotalObjects = n;}
    /* File IO */
    void record_frame();
    void write_frame_header();
    /* World controls */
    void createObjects(unsigned n) {for (unsigned i=0; i<n; i++) createObject();}
    void createObject();
    void advance_full_step();
    void advance_half_step();
    void exchange_ghost_objects();
    void exchange_objects();
    void detect_collisions();

    void set_neighbors(std::vector<Slice_serial *> slices);
    void clear_ghost_objects();
    void addGhostObject(Object * obj) {ghostObjects.push_back(obj);}
    void addObject(Object * obj) {objects.push_back(obj);}
private:
    unsigned rank;
    unsigned nSlices;
    unsigned nTotalObjects; // count of all objects in the world
    std::ofstream * fileHandle;
    double overlapWidth;
    Bounds bounds;
    Bounds innerBounds;
    Neighbors_rank neighbors_rank;
    Neighbors_pointer neighbors;
    std::vector<Object *> objects;
    std::vector<Object *> ghostObjects;

    /* object IDs */
    unsigned nextObjectID;
    unsigned minID;
    unsigned maxID; 
    unsigned nextFrameID;

    /* Private methods */

    void detect_collisions_world_boundaries();
    void handle_collision(Object * obj1, Object * obj2);


};

