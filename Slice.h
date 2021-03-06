#pragma once

#include "vec3.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <mpi.h>

#include "Object.h"

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

class Slice {
public:
    Slice(int r, int s, MPI_File * fh, double worldSize[3], unsigned sliceRows[3]);
    ~Slice();
    void setTotalObjects(unsigned n) {nTotalObjects = n;}
    /* File IO */
    void record_frame();
    void write_frame_header();
    /* World controls */
    void createObjects(int n) {for (int i=0; i<n; i++) createObject();}
    void createObject();
    void advance_full_step();
    void advance_half_step();
    void exchange_ghost_objects();
    void exchange_objects();
    void detect_collisions();
private:
    int rank;
    int nSlices;
    unsigned nTotalObjects; // count of all objects in the world
    MPI_File * fileHandle;
    double overlapWidth;
    Bounds bounds;
    Bounds innerBounds;
    Neighbors neighbors;
    std::vector<Object *> objects;
    std::vector<Object *> ghostObjects;

    /* object IDs */
    unsigned nextObjectID;
    int minID;
    int maxID; 
    unsigned nextFrameID;

    /* Private methods */
    void store_object_state(std::vector<unsigned> * idBuffer, std::vector<double> * stateBuffer, Object * o);
    void send_objects(std::vector<unsigned> * idBuffer, std::vector<double> * stateBuffer, unsigned targetRank, MPI_Request * requestList);
    void receive_objects(unsigned sourceRank, std::vector<Object *> * objList);

    void detect_collisions_world_boundaries();
    void handle_collision(Object * obj1, Object * obj2);


};

