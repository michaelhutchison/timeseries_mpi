#include "Slice.h"

/* Constructor */
Slice::Slice(int r, int s, MPI_File * fh, double worldSize[3]) {
    rank = r;
    nSlices = s;
    fileHandle = fh;

    // initialize world boundaries
    // -- This might be better to initialize outside this class.
    Bounds worldBounds;
    worldBounds.nx = -worldSize[0] * 0.5;
    worldBounds.px =  worldSize[0] * 0.5;
    worldBounds.ny = -worldSize[1] * 0.5;
    worldBounds.py =  worldSize[1] * 0.5;
    worldBounds.nz = -worldSize[2] * 0.5;
    worldBounds.pz =  worldSize[2] * 0.5;

    // initialize slice boundaries 
    // -- Each slice has a portion of the x-range 
    // -- Each slice covers the entire y-range and z-range
    bounds = worldBounds;
    double sliceWidth = (worldBounds.px - worldBounds.nx)/nSlices;
    bounds.nx = worldBounds.nx + (rank * sliceWidth);
    bounds.px = bounds.nx + sliceWidth;
    overlapWidth = 0.3;
    innerBounds = bounds;
    if (rank > 0) innerBounds.nx += overlapWidth;
    if (rank < nSlices-1) innerBounds.px -= overlapWidth;
    
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
    nTotalObjects = 0;
}
Slice::~Slice() {
    // Release dynamically allocated memory
    for (int i=0; i<objects.size(); i++) 
        delete objects[i];
}
void Slice::write_frame_header() {
    //std::cout << "I got " << nTotalObjects << " total objects. I have " << objects.size() << std::endl;
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
void Slice::record_frame() {
    // Save the current frame to the animation record

    // Collect the number of objects. 
    // This block is commented out because the number of objects is static.
    /*
    unsigned nLocalObjects = objects.size();
    MPI_Reduce( &nLocalObjects,
                &nTotalObjects,
                1,
                MPI_UNSIGNED,
                MPI_SUM,
                0,
                MPI_COMM_WORLD);
    */
    // process with rank 0 writes the frame header
    if (rank == 0) {
        write_frame_header();
    }
    // Let all processes catch up
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_seek_shared(*fileHandle, 0, MPI_SEEK_END);
    
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
    double vMax = 0.5;
    o->setVelocity( RANDOMDOUBLE(-vMax, vMax),
                    RANDOMDOUBLE(-vMax, vMax),
                    RANDOMDOUBLE(-vMax, vMax) );

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
void Slice::send_objects(std::vector<unsigned> * idBuffer, std::vector<double> * stateBuffer, unsigned targetRank, MPI_Request * requestList) {
    unsigned doublesPerObject = 10;
    // Send number of objects
    unsigned nObjects = idBuffer->size();
    MPI_Isend(&nObjects,
             1,
             MPI_UNSIGNED,
             targetRank,
             1,
             MPI_COMM_WORLD,
             &requestList[0]);
    // Send object ids
    std::vector<unsigned> &ids = *idBuffer;
    MPI_Isend(&ids[0], 
             nObjects,
             MPI_UNSIGNED,
             targetRank,
             2,
             MPI_COMM_WORLD,
             &requestList[1]);
    // Send object states
    std::vector<double> &states = *stateBuffer;
    MPI_Isend(&states[0], 
             nObjects * doublesPerObject,
             MPI_DOUBLE,
             targetRank,
             3,
             MPI_COMM_WORLD,
             &requestList[2]);
}

/* receive_objects
 *   Accepts object ids and states for objects which are
 *   being transferred from another slice.
 *   Creates and ingests new objects once data has been received.
 */
void Slice::receive_objects(unsigned sourceRank, std::vector<Object_mpi *> * objList) {
    unsigned doublesPerObject = 10;
    // Receive number of objects
    unsigned nObjects;
    MPI_Request request;
    MPI_Irecv(&nObjects,
             1,
             MPI_UNSIGNED,
             sourceRank,
             1,
             MPI_COMM_WORLD,
             &request);
    MPI_Wait(&request, MPI_STATUSES_IGNORE);
    // Receive object ids
    unsigned * idBuffer = new unsigned[nObjects];
    MPI_Irecv(idBuffer,
             nObjects,
             MPI_UNSIGNED,
             sourceRank,
             2,
             MPI_COMM_WORLD,
             &request);
    MPI_Wait(&request, MPI_STATUSES_IGNORE);
    // Receive object states
    double * stateBuffer = new double[nObjects * doublesPerObject];
    MPI_Irecv(stateBuffer,
             nObjects * doublesPerObject,
             MPI_DOUBLE,
             sourceRank,
             3,
             MPI_COMM_WORLD,
             &request);
    MPI_Wait(&request, MPI_STATUSES_IGNORE);
    // Ingest new objects
    for (int i=0; i<nObjects; i++) {
        Object_mpi * o = new Object_mpi;
        objList->push_back(o);
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
void Slice::exchange_ghost_objects() {
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

    // Results objects
    MPI_Request leftSendRequestList[3];
    MPI_Request rightSendRequestList[3];

    // Clear the ghost object list
    ghostObjects.clear();

    // Send to negative-x neighbor
    if (neighbors.nx != -1) {
        idBuffer.clear();
        stateBuffer.clear();
        // Produce buffer of objects to transfer to nx neighbor
        for (std::vector<Object_mpi *>::iterator itr = objects.begin(); itr != objects.end(); ++itr) {
            if ((*itr)->x() < innerBounds.nx) {
                // Copy the object's id and state into transfer buffers
                store_object_state(&idBuffer, &stateBuffer, (*itr));
            }
        }
        // Transfer the set of objects to the left (negative-x) neighbor
        send_objects(&idBuffer, &stateBuffer, rank-1, leftSendRequestList);
    }
    // Send to positive-x neighbor
    if (neighbors.px != -1) {
        idBuffer.clear();
        stateBuffer.clear();
        // Produce buffer of objects to transfer to px neighbor
        for (std::vector<Object_mpi *>::iterator itr = objects.begin(); itr != objects.end(); ++itr) {
            if ((*itr)->x() >= innerBounds.px) {
                // Copy the object's id and state into transfer buffers
                store_object_state(&idBuffer, &stateBuffer, (*itr));
            }
        }
        // Transfer the set of objects to the right (positive-x) neighbor
        send_objects(&idBuffer, &stateBuffer, rank+1, rightSendRequestList);
        
    }

   // Receive from negative-x neighbors
    if (neighbors.nx != -1) {
        receive_objects(rank-1, &ghostObjects);
    }
    // Receive from positive-x neighbors
    if (neighbors.px != -1) {
        receive_objects(rank+1, &ghostObjects);
    }
    // Wait for communication to be complete
    if (neighbors.nx != -1) { MPI_Waitall(3, leftSendRequestList, MPI_STATUSES_IGNORE); }
    if (neighbors.px != -1) { MPI_Waitall(3, rightSendRequestList, MPI_STATUSES_IGNORE); }
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

    // Results objects
    MPI_Request leftSendRequestList[3];
    MPI_Request rightSendRequestList[3];

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
        send_objects(&idBuffer, &stateBuffer, rank-1, leftSendRequestList);
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
        send_objects(&idBuffer, &stateBuffer, rank+1, rightSendRequestList);
    }

   // Receive from negative-x neighbors
    if (neighbors.nx != -1) {
        receive_objects(rank-1, &objects);
    }
     // Receive from positive-x neighbors
    if (neighbors.px != -1) {
        receive_objects(rank+1, &objects);
    }
 
    // Wait for communication to be complete
    if (neighbors.nx != -1) { MPI_Waitall(3, leftSendRequestList, MPI_STATUSES_IGNORE); }
    if (neighbors.px != -1) { MPI_Waitall(3, rightSendRequestList, MPI_STATUSES_IGNORE); }
 
}

void Slice::handle_collision(Object_mpi * obj1, Object_mpi * obj2) {

    // Get position vectors
    Vec3 obj1Position = obj1->getPosition();
    Vec3 obj2Position = obj2->getPosition();
    // Get velocity vectors
    Vec3 obj1Velocity = obj1->getVelocity();
    Vec3 obj2Velocity = obj2->getVelocity();
    // Save speed magnitude
    double speed1 = obj1Velocity.length();
    double speed2 = obj2Velocity.length();
    // calculate normal center-to-center vector and magnitude of overlap
    Vec3 sphereNormal = obj2Position - obj1Position;
    double overlap = (obj1->getRadius() + obj2->getRadius()) - sphereNormal.length();
    sphereNormal.normalize();
    // calculate object 1 new velocity 
    double dotProd = obj1Velocity.dot(sphereNormal);
    Vec3 obj1NewVelocity = obj1Velocity - 2*dotProd*sphereNormal;
    // calculate object 2 new velocity 
    dotProd = obj2Velocity.dot(-1*sphereNormal);
    Vec3 obj2NewVelocity = obj2Velocity + 2*dotProd*sphereNormal;
    // Fix speeds
    obj1NewVelocity *= speed2/obj1NewVelocity.length();
    obj2NewVelocity *= speed1/obj2NewVelocity.length();
   
    obj1->setVelocity(obj1NewVelocity);
    obj2->setVelocity(obj2NewVelocity);
    // Fix positions so spheres no longer overlap
    obj1->delta_position(sphereNormal * overlap * -0.55);
    obj2->delta_position(sphereNormal * overlap * 0.55);


}
void Slice::detect_collisions() {
    // Handle collisions with ghost objects
    for (std::vector<Object_mpi *>::iterator itr = objects.begin(); itr != objects.end(); ++itr) {
        for (std::vector<Object_mpi *>::iterator ghostItr = ghostObjects.begin(); ghostItr != ghostObjects.end(); ++ghostItr) {
            // Get distance between centers of spheres
            double distance = sqrt( pow((*itr)->x()-(*ghostItr)->x(), 2) + 
                                    pow((*itr)->y()-(*ghostItr)->y(), 2) +
                                    pow((*itr)->z()-(*ghostItr)->z(), 2) );
            if (distance < 2*(*itr)->getRadius()) {
                handle_collision((*itr), (*ghostItr));
                //std::cout << "BOUNCE" << std::endl;
            }
        }    
    }
    // Handle collisions with local objects
    //double minDistance = 20.0;
    for (std::vector<Object_mpi *>::iterator itr = objects.begin(); itr != objects.end(); ++itr) {
        for (std::vector<Object_mpi *>::iterator otherItr = itr; otherItr != objects.end(); ++otherItr) {
            if ((*itr)->getID() != (*otherItr)->getID()) {
                // Get distance between centers of spheres
                double distance = sqrt( pow((*itr)->x()-(*otherItr)->x(), 2) + 
                                        pow((*itr)->y()-(*otherItr)->y(), 2) +
                                        pow((*itr)->z()-(*otherItr)->z(), 2) );
                //if (distance < minDistance) minDistance = distance;
                //if (rank == 0) std::cout << "radius: " << (*itr)->getRadius() << std::endl;
                if (distance < 2*(*itr)->getRadius()) {
                    handle_collision((*itr), (*otherItr));
                    //std::cout << "LOCAL BOUNCE" << std::endl;
                }
            }
        }    
    }
    //if (rank == 0) std::cout << "min distance: " << minDistance << std::endl;
    detect_collisions_world_boundaries();
}
void Slice::detect_collisions_world_boundaries() {
    
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
