#include "Slice.h"

/* Constructor */
Slice::Slice(int r, int s, MPI_File * fh) {
    rank = r;
    nSlices = s;
    fileHandle = fh;

    // initialize world boundaries
    // -- This might be better to initialize outside this class.
    Bounds worldBounds;
    worldBounds.nx = -6.0;
    worldBounds.px =  6.0;
    worldBounds.ny = -6.0;
    worldBounds.py =  6.0;
    worldBounds.nz = -6.0;
    worldBounds.pz =  6.0;

    // initialize slice boundaries 
    // -- Each slice has a portion of the x-range 
    // -- Each slice covers the entire y-range and z-range
    bounds = worldBounds;
    double sliceWidth = (worldBounds.px - worldBounds.nx)/nSlices;
    bounds.nx = worldBounds.nx + (rank * sliceWidth);
    bounds.px = bounds.nx + sliceWidth;
    
    // Initialize neighbors
    // -- value of -1 indicates no neighbor at that boundary
    neighbors.nx = -1;
    neighbors.px = -1;
    neighbors.ny = -1;
    neighbors.py = -1;
    neighbors.nz = -1;
    neighbors.pz = -1;
    if (rank > 0) neighbors.nx = rank-1;
    if (rank < nSlices-1) neighbors.px = rank+1;

    // Calculate id range
    // Each slice assigns an id to the objects it creates.
    // -- upper 4 bits are reserved for slice id -- MAXIMUM OF 16 SLICES
    // -- lower 28 bits are reserved for object id -- MAXIMUM OF 268,435,455 OBJECTS PER SLICE
    // -- ID range: [rank*2^28] to [((rank+1)*2^28)-1]
    minID = rank << 28;
    maxID = ((rank+1) << 28)-1;
    nextObjectID = minID;
    nextFrameID = 0;

}
Slice::~Slice() {
    // Release dynamically allocated memory
    for (int i=0; i<objects.size(); i++) 
        delete objects[i];
}
void Slice::createObject() {
    // Do nothing if this slice has used its ID range
    if (nextObjectID > maxID) return;
    // Create new object, assign ID
    Object_mpi * o = new Object_mpi;
    objects.push_back(o);
    o->setID(nextObjectID);
    nextObjectID++;
    // Set object's state
    o->setPosition(  RANDOMDOUBLE(bounds.nx, bounds.px),
                     RANDOMDOUBLE(bounds.ny, bounds.py),
                     RANDOMDOUBLE(bounds.nz, bounds.pz) );
    o->setRotationVector( RANDOMDOUBLE(-1.0, 1.0),
                          RANDOMDOUBLE(-1.0, 1.0),
                          RANDOMDOUBLE(-1.0, 1.0) );
    o->setRotationAngle( RANDOMDOUBLE(0, 2*PI));
    o->setRotationVelocity( RANDOMDOUBLE(-5,5));
    double vMax = 0.1;
    o->setVelocity( RANDOMDOUBLE(-vMax, vMax),
                    RANDOMDOUBLE(-vMax, vMax),
                    RANDOMDOUBLE(-vMax, vMax) );

    if (nextObjectID == 1) {
        std::cout << "Object id: " << o->getID() << std::endl;
        std::cout << "object position: " << o->x() << "  " << o->y() << "  " << o->z() << std::endl;
    }
}

void Slice::record_frame() {
    // Save the current frame to the animation record
    //timeseries.record_frame(objects);

    // Collect the number of objects. 
    unsigned nTotalObjects;
    unsigned nLocalObjects = objects.size();
    MPI_Reduce( &nLocalObjects,
                &nTotalObjects,
                1,
                MPI_UNSIGNED,
                MPI_SUM,
                0,
                MPI_COMM_WORLD);
    
    // process with rank 0 writes the frame header
    if (rank == 0) {
        std::cout << "I got " << nTotalObjects << " total objects. I have " << objects.size() << std::endl;
        unsigned short objectLength = Object_mpi::objectLengthInBytes;
        unsigned long frameSizeInBytes = nTotalObjects * objectLength;
        MPI_Status status;
        MPI_File_seek(*fileHandle,
                         0,
                         MPI_SEEK_END);
 
        MPI_File_write(*fileHandle,
                       &frameSizeInBytes,
                       1,
                       MPI_UNSIGNED_LONG,
                       &status);
        MPI_File_write(*fileHandle,
                       &nextFrameID,
                       1,
                       MPI_UNSIGNED,
                       &status);
        MPI_File_write(*fileHandle,
                       &nTotalObjects,
                       1,
                       MPI_UNSIGNED,
                       &status);
        MPI_File_write(*fileHandle,
                       &objectLength,
                       1,
                       MPI_UNSIGNED_SHORT,
                       &status);

        nextFrameID++;
    }
    // Let all processes catch up
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_seek_shared(*fileHandle,
                         0,
                         MPI_SEEK_END);
    
    // Create buffer containing all object states
    unsigned totalBytes = objects.size() * Object_mpi::objectLengthInBytes;
    char * objectBuffer = new char[totalBytes];
    memset(objectBuffer, 0, sizeof(objectBuffer));
    for (int i=0; i<objects.size(); i++) {
        objects[i]->get_bytes(&objectBuffer[i*60]);
    }
    // Write the buffer
    MPI_Status status;
    MPI_File_write_ordered(*fileHandle,
                        objectBuffer,
                        totalBytes,
                        MPI_BYTE,
                        &status);
}

void Slice::advance_full_step() {
    /* Advances position and velocity a half-step for each object */
    for (int i=0; i<objects.size(); i++) {
        objects[i]->advance_full_step();
    }   
}
void Slice::advance_half_step() {
    /* Advances position and velocity a half-step for each object */
    for (int i=0; i<objects.size(); i++) {
        objects[i]->advance_half_step();
    }   
}

/* store_object_state
 *   copies an object's id and crucial state information to buffers
 *   in preparation to transfer the object to another slice.
 */
void Slice::store_object_state(std::vector<unsigned> * idBuffer, std::vector<double> * stateBuffer, Object_mpi * o) {
    // Save object id
    idBuffer->push_back(o->getID());
    // Save object state
    Vec3 position = o->getPosition();
    Vec3 velocity = o->getVelocity();
    Vec3 rotationVector = o->getRotationVector();
    stateBuffer->push_back(position.x);
    stateBuffer->push_back(position.y);
    stateBuffer->push_back(position.z);
    stateBuffer->push_back(velocity.x);
    stateBuffer->push_back(velocity.y);
    stateBuffer->push_back(velocity.z);
    stateBuffer->push_back(rotationVector.x);
    stateBuffer->push_back(rotationVector.y);
    stateBuffer->push_back(rotationVector.z);
    stateBuffer->push_back(o->getRotationAngle());
}

/* send_objects
 *   Transfers object ids and states to another slice.
 *   These objects should already have been removed from the 
 *   list of objects owned by the sending slice.
 */
void Slice::send_objects(std::vector<unsigned> * idBuffer, std::vector<double> * stateBuffer, unsigned targetRank) {
    int doublesPerObject = 10;
    // Send number of objects
    unsigned nObjects = idBuffer->size();
    MPI_Send(&nObjects,
             1,
             MPI_UNSIGNED,
             targetRank,
             0,
             MPI_COMM_WORLD);
    if (nObjects == 0) return;
    // Send object ids
    std::vector<unsigned> &ids = *idBuffer;
    MPI_Send(&ids[0], 
             idBuffer->size(),
             MPI_UNSIGNED,
             targetRank,
             0,
             MPI_COMM_WORLD);
    // Send object states
    std::vector<double> &states = *stateBuffer;
    MPI_Send(&states[0], 
             idBuffer->size() * doublesPerObject,
             MPI_DOUBLE,
             targetRank,
             0,
             MPI_COMM_WORLD);
}

/* receive_objects
 *   Accepts object ids and states for objects which are
 *   being transferred from another slice.
 *   Creates and ingests new objects once data has been received.
 */
void Slice::receive_objects(unsigned sourceRank) {
    MPI_Status status;
    int doublesPerObject = 10;
    // Receive number of objects
    unsigned nObjects;
    MPI_Recv(&nObjects,
             1,
             MPI_UNSIGNED,
             sourceRank,
             0,
             MPI_COMM_WORLD,
             &status);
    if (nObjects == 0) return;
    // Receive object ids
    unsigned * idBuffer = new unsigned[nObjects];
    MPI_Recv(idBuffer,
             nObjects,
             MPI_UNSIGNED,
             sourceRank,
             0,
             MPI_COMM_WORLD,
             &status);
    // Receive object states
    double * stateBuffer = new double[nObjects * doublesPerObject];
    MPI_Recv(stateBuffer,
             nObjects * doublesPerObject,
             MPI_DOUBLE,
             sourceRank,
             0,
             MPI_COMM_WORLD,
             &status);
    // Ingest new objects
    for (int i=0; i<nObjects; i++) {
        Object_mpi * o = new Object_mpi;
        objects.push_back(o);
        o->setID(idBuffer[i]);
        Vec3 position(stateBuffer[i*doublesPerObject + 0],
                      stateBuffer[i*doublesPerObject + 1],
                      stateBuffer[i*doublesPerObject + 2]);
        Vec3 velocity(stateBuffer[i*doublesPerObject + 3],
                      stateBuffer[i*doublesPerObject + 4],
                      stateBuffer[i*doublesPerObject + 5]);
        Vec3 rotation(stateBuffer[i*doublesPerObject + 6],
                      stateBuffer[i*doublesPerObject + 7],
                      stateBuffer[i*doublesPerObject + 8]);
        double angle = stateBuffer[i*doublesPerObject + 9];
        o->setPosition(position);
        o->setVelocity(velocity);
        o->setRotationVector(rotation);
        o->setRotationAngle(angle);
    }
}

void Slice::exchange_objects() {
    // Any objects which have crossed a boundary with another slice
    // are transferred to that slice. 
    // Any objects which have crossed into this slice are
    // recieved from their previous slice.
    // NOTE: Any objects whose anchor points are exactly on a slice's 
    //       negative-x boundary are owned by that slice. 
    //       An object on a slice's positive-x boundary are owned 
    //       by their neighber with rank+1.
    //
    // For each object, two buffers are sent:
    // 1) object ids
    //   unsigned            object ID 
    // 2) object states
    //   double              x position
    //   double              y position
    //   double              z position
    //   double              x velocity
    //   double              y velocity
    //   double              z velocity
    //   double              x of rotation vector
    //   double              y of rotation vector
    //   double              z of rotation vector
    //   double              angle of rotation
    
    // Create buffer of objects to transfer
    std::vector<unsigned> idBuffer;
    std::vector<double> stateBuffer;
    int doublesPerObject = 10;

    // Send to negative-x neighbor
    if (neighbors.nx != -1) {
        idBuffer.clear();
        stateBuffer.clear();
        // Produce buffer of objects to transfer to nx neighbor
        for (std::vector<Object_mpi *>::iterator itr = objects.begin(); itr != objects.end(); ) {
            if ((*itr)->x() < bounds.nx) {
                // Copy the object's id and state into transfer buffers
                store_object_state(&idBuffer, &stateBuffer, (*itr));
                // Delete the object reference from the local list
                itr = objects.erase(itr);
            } else {
                // Advance the iterator
                ++itr;
            }
        }
        // Transfer the set of objects to the neighbor
        send_objects(&idBuffer, &stateBuffer, rank-1);
    }
    // Receive from positive-x neighbors
    if (neighbors.px != -1) {
        receive_objects(rank+1);
    }
    // Send to positive-x neighbor
    if (neighbors.px != -1) {
        idBuffer.clear();
        stateBuffer.clear();
        // Produce buffer of objects to transfer to nx neighbor
        for (std::vector<Object_mpi *>::iterator itr = objects.begin(); itr != objects.end(); ) {
            if ((*itr)->x() >= bounds.px) {
                // Copy the object's id and state into transfer buffers
                store_object_state(&idBuffer, &stateBuffer, (*itr));
                // Delete the object reference from the local list
                itr = objects.erase(itr);
            } else {
                // Advance the iterator
                ++itr;
            }
        }
        // Transfer the set of objects to the neighbor
        send_objects(&idBuffer, &stateBuffer, rank+1);
    }
    // Receive from negative-x neighbors
    if (neighbors.nx != -1) {
        receive_objects(rank-1);
    }
 


}
void Slice::handle_collisions() {
    // NOTE: only handles collisions with world boundaries.
    for (std::vector<Object_mpi *>::iterator itr = objects.begin(); itr != objects.end(); ++itr) {
        Vec3 dp;
        Vec3 dv(1,1,1);

        // x boundaries
        if (neighbors.px == -1 && (*itr)->x() + (*itr)->getRadius() > bounds.px) {
            dp.x = bounds.px - ((*itr)->x() + (*itr)->getRadius());
            dv.x = -1;
        } 
        else if (neighbors.nx == -1 && (*itr)->x() - (*itr)->getRadius() < bounds.nx) {
            dp.x = bounds.nx - ((*itr)->x() - (*itr)->getRadius()); 
            dv.x = -1;
        }
        // y boundaries
        if ((*itr)->y() + (*itr)->getRadius() > bounds.py) {
            dp.y = bounds.py - ((*itr)->y() + (*itr)->getRadius());
            dv.y = -1;
        } 
        else if ((*itr)->y() - (*itr)->getRadius() < bounds.ny) {
            dp.y = bounds.ny - ((*itr)->y() - (*itr)->getRadius()); 
            dv.y = -1;
        }
        // z boundaries
        if ((*itr)->z() + (*itr)->getRadius() > bounds.pz) {
            dp.z = bounds.pz - ((*itr)->z() + (*itr)->getRadius());
            dv.z = -1; 
        } 
        else if ((*itr)->z() - (*itr)->getRadius() < bounds.nz) {
            dp.z = bounds.nz - ((*itr)->z() - (*itr)->getRadius()); 
            dv.z = -1;
        }
        // Update object position
        (*itr)->delta_position(dp);
        (*itr)->multiply_velocity(dv);
    }
}