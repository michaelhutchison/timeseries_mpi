#pragma once

#include "vec3.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <stdlib.h>
#include <time.h>
#include <cmath>

#include "Object_mpi.h"

#define RANDOMDOUBLE(min,max) (  ((double)rand() / RAND_MAX) * (max-min) + min  ) 
#ifndef PI
#define PI 3.14159265
#endif

struct Bounds {
    double nx; // negative x boundary
    double px; // positive x boundary
    double ny; // negative y boundary
    double py; // positive y boundary
    double nz; // negative z boundary
    double pz; // positive z boundary
};
struct Neighbors {
    // Holds the rank of the neighbor at each boundary.
    // Value of -1 indicates no neighbor on that boundary.
    int nx; // negative x neighbor
    int px; // positive x neighbor
    int ny; // negative y neighbor
    int py; // positive y neighbor
    int nz; // negative z neighbor
    int pz; // positive z neighbor
};

class World {
public:
    World(std::ofstream * fh, double worldSize[3]);
    ~World();
    /* File IO */
    void record_frame();
    void write_frame_header();
    /* World controls */
    void createObjects(int n) {for (int i=0; i<n; i++) createObject();}
    void createObject();
    void advance_full_step();
    void advance_half_step();
    void detect_collisions();
private:
    std::ofstream * fileHandle;
    Bounds bounds;
    std::vector<Object_mpi *> objects;
    /* boundary cases */

    /* object IDs */
    unsigned nextObjectID;
    unsigned minID;
    unsigned maxID; 
    unsigned nextFrameID;

    /* Private methods */
    void detect_collisions_world_boundaries();
    void handle_collision(Object_mpi * obj1, Object_mpi * obj2);


};

