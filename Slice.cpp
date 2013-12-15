#include "Slice.h"

/* Constructor */
Slice::Slice(int r, int s, MPI_File * fh, double worldSize[3], unsigned sliceRows[3]) {
    rank = r;
    nSlices = s;
    fileHandle = fh;

    unsigned slicesX = sliceRows[0];
    unsigned slicesY = sliceRows[1];
    unsigned slicesZ = sliceRows[2];

    // initialize world boundaries
    // -- This might be better to initialize outside this class.
    Bounds worldBounds;
    worldBounds.nx = -worldSize[0] * 0.5;
    worldBounds.px =  worldSize[0] * 0.5;
    worldBounds.ny = -worldSize[1] * 0.5;
    worldBounds.py =  worldSize[1] * 0.5;
    worldBounds.nz = -worldSize[2] * 0.5;
    worldBounds.pz =  worldSize[2] * 0.5;

    // determine position in rows and columns based on rank
    unsigned rowX = rank % slicesX;
    unsigned rowY = (rank/slicesX) % slicesY;
    unsigned rowZ = (rank/(slicesX*slicesY)) % slicesZ;
    unsigned minInRowX = (rowY * slicesX) + (rowZ*slicesX*slicesY);
    unsigned maxInRowX = minInRowX + slicesX - 1;
    unsigned minInRowY = rowX % slicesX + (rowZ*slicesX*slicesY);
    unsigned maxInRowY = minInRowY + (slicesX*slicesY) - slicesX;
    unsigned minInRowZ = (rowY*slicesX) + rowX;
    unsigned maxInRowZ = minInRowZ + (slicesX*slicesY*(slicesZ-1));

    // initialize slice boundaries 
    // -- Each slice has a portion of the x-range 
    // -- Each slice covers the entire y-range and z-range
    bounds = worldBounds;
    double sliceWidthX = (worldBounds.px - worldBounds.nx)/slicesX;
    double sliceWidthY = (worldBounds.py - worldBounds.ny)/slicesY;
    double sliceWidthZ = (worldBounds.pz - worldBounds.nz)/slicesZ;
    bounds.nx = worldBounds.nx + (rowX * sliceWidthX);
    bounds.px = bounds.nx + sliceWidthX;
    bounds.ny = worldBounds.ny + (rowY * sliceWidthY);
    bounds.py = bounds.ny + sliceWidthY;
    bounds.nz = worldBounds.nz + (rowZ * sliceWidthZ);
    bounds.pz = bounds.nz + sliceWidthZ;
    // inner boundaries
    overlapWidth = 0.3;
    innerBounds = bounds;
    if (rank > minInRowX) innerBounds.nx += overlapWidth;
    if (rank < maxInRowX) innerBounds.px -= overlapWidth;
    if (rank > minInRowY) innerBounds.ny += overlapWidth;
    if (rank < maxInRowY) innerBounds.py -= overlapWidth;
    if (rank > minInRowZ) innerBounds.nz -= overlapWidth;
    if (rank < maxInRowZ) innerBounds.pz -= overlapWidth;

    // Initialize neighbors
    // -- value of -1 indicates no neighbor at that boundary
    neighbors.nx = -1;
    neighbors.px = -1;
    neighbors.ny = -1;
    neighbors.py = -1;
    neighbors.nz = -1;
    neighbors.pz = -1;
    if (rank > minInRowX) neighbors.nx = rank-1;
    if (rank < maxInRowX) neighbors.px = rank+1;
    if (rank > minInRowY) neighbors.ny = rank-slicesX;
    if (rank < maxInRowY) neighbors.py = rank+slicesX;
    if (rank > minInRowZ) neighbors.nz = rank-(slicesX*slicesY);
    if (rank < maxInRowZ) neighbors.pz = rank+(slicesX*slicesY);


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
    unsigned short objectLength = Object::objectLengthInBytes;
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
    unsigned totalBytes = objects.size() * Object::objectLengthInBytes;
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
    Object * o = new Object;
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
void Slice::store_object_state(std::vector<unsigned> * idBuffer, std::vector<double> * stateBuffer, Object * o) {
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
void Slice::receive_objects(unsigned sourceRank, std::vector<Object *> * objList) {
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
        Object * o = new Object;
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
    std::vector<unsigned> idBufferNX;
    std::vector<unsigned> idBufferPX;
    std::vector<unsigned> idBufferNY;
    std::vector<unsigned> idBufferPY;
    std::vector<unsigned> idBufferNZ;
    std::vector<unsigned> idBufferPZ;
    std::vector<double> stateBufferNX;
    std::vector<double> stateBufferPX;
    std::vector<double> stateBufferNY;
    std::vector<double> stateBufferPY;
    std::vector<double> stateBufferNZ;
    std::vector<double> stateBufferPZ;
    int doublesPerObject = 10;

    // Results objects
    MPI_Request nxSendRequestList[3];
    MPI_Request pxSendRequestList[3];
    MPI_Request nySendRequestList[3];
    MPI_Request pySendRequestList[3];
    MPI_Request nzSendRequestList[3];
    MPI_Request pzSendRequestList[3];

    // Clear the ghost object list
    for (unsigned i=0; i<ghostObjects.size(); i++)
        delete ghostObjects[i];
    ghostObjects.clear();

    // Examine each object
    for (std::vector<Object *>::iterator itr = objects.begin(); itr != objects.end(); ++itr) {
        if (neighbors.nx != -1 && (*itr)->x() < innerBounds.nx) {
            // Copy the object's id and state into transfer buffers
            store_object_state(&idBufferNX, &stateBufferNX, (*itr));
        } 
        if (neighbors.px != -1 && (*itr)->x() > innerBounds.px) {
            // Copy the object's id and state into transfer buffers
            store_object_state(&idBufferPX, &stateBufferPX, (*itr));
        }
        if (neighbors.ny != -1 && (*itr)->y() < innerBounds.ny) {
            // Copy the object's id and state into transfer buffers
            store_object_state(&idBufferNY, &stateBufferNY, (*itr));
        }
        if (neighbors.py != -1 && (*itr)->y() > innerBounds.py) {
            // Copy the object's id and state into transfer buffers
            store_object_state(&idBufferPY, &stateBufferPY, (*itr));
        } 
        if (neighbors.nz != -1 && (*itr)->z() < innerBounds.nz) {
            // Copy the object's id and state into transfer buffers
            store_object_state(&idBufferNZ, &stateBufferNZ, (*itr));
        }
        if (neighbors.pz != -1 && (*itr)->z() > innerBounds.pz) {
            // Copy the object's id and state into transfer buffers
            store_object_state(&idBufferPZ, &stateBufferPZ, (*itr));
        } 
 
    }

    // Send object lists
    if (neighbors.nx != -1) send_objects(&idBufferNX, &stateBufferNX, neighbors.nx, nxSendRequestList);
    if (neighbors.px != -1) send_objects(&idBufferPX, &stateBufferPX, neighbors.px, pxSendRequestList);
    if (neighbors.ny != -1) send_objects(&idBufferNY, &stateBufferNY, neighbors.ny, nySendRequestList);
    if (neighbors.py != -1) send_objects(&idBufferPY, &stateBufferPY, neighbors.py, pySendRequestList);
    if (neighbors.nz != -1) send_objects(&idBufferNZ, &stateBufferNZ, neighbors.nz, nzSendRequestList);
    if (neighbors.pz != -1) send_objects(&idBufferPZ, &stateBufferPZ, neighbors.pz, pzSendRequestList);

    // Receive from negative-x neighbors
    if (neighbors.nx != -1) {
        receive_objects(rank-1, &ghostObjects);
    }
    // Receive from positive-x neighbors
    if (neighbors.px != -1) {
        receive_objects(rank+1, &ghostObjects);
    }
    // Receive from negative-y neighbors
    if (neighbors.ny != -1) {
        receive_objects(neighbors.ny, &ghostObjects);
    }
    // Receive from positive-y neighbors
    if (neighbors.py != -1) {
        receive_objects(neighbors.py, &ghostObjects);
    }
    // Receive from negative-z neighbors
    if (neighbors.nz != -1) {
        receive_objects(neighbors.nz, &ghostObjects);
    }
    // Receive from positive-z neighbors
    if (neighbors.pz != -1) {
        receive_objects(neighbors.pz, &ghostObjects);
    }
 
    // Wait for communication to be complete
    if (neighbors.nx != -1) { MPI_Waitall(3, nxSendRequestList, MPI_STATUSES_IGNORE); }
    if (neighbors.px != -1) { MPI_Waitall(3, pxSendRequestList, MPI_STATUSES_IGNORE); }
    if (neighbors.ny != -1) { MPI_Waitall(3, nySendRequestList, MPI_STATUSES_IGNORE); }
    if (neighbors.py != -1) { MPI_Waitall(3, pySendRequestList, MPI_STATUSES_IGNORE); }
    if (neighbors.nz != -1) { MPI_Waitall(3, nzSendRequestList, MPI_STATUSES_IGNORE); }
    if (neighbors.pz != -1) { MPI_Waitall(3, pzSendRequestList, MPI_STATUSES_IGNORE); }
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
    std::vector<unsigned> idBufferNX;
    std::vector<unsigned> idBufferPX;
    std::vector<unsigned> idBufferNY;
    std::vector<unsigned> idBufferPY;
    std::vector<unsigned> idBufferNZ;
    std::vector<unsigned> idBufferPZ;
    std::vector<double> stateBufferNX;
    std::vector<double> stateBufferPX;
    std::vector<double> stateBufferNY;
    std::vector<double> stateBufferPY;
    std::vector<double> stateBufferNZ;
    std::vector<double> stateBufferPZ;
    int doublesPerObject = 10;

    // Results objects
    MPI_Request nxSendRequestList[3];
    MPI_Request pxSendRequestList[3];
    MPI_Request nySendRequestList[3];
    MPI_Request pySendRequestList[3];
    MPI_Request nzSendRequestList[3];
    MPI_Request pzSendRequestList[3];

    // Examine each object
    for (std::vector<Object *>::iterator itr = objects.begin(); itr != objects.end(); ) {
        if (neighbors.nx != -1 && (*itr)->x() < bounds.nx) {
            // Copy the object's id and state into transfer buffers
            store_object_state(&idBufferNX, &stateBufferNX, (*itr));
            // Delete the object reference from the local list
            itr = objects.erase(itr);
        } else if (neighbors.px != -1 && (*itr)->x() > bounds.px) {
            // Copy the object's id and state into transfer buffers
            store_object_state(&idBufferPX, &stateBufferPX, (*itr));
            // Delete the object reference from the local list
            itr = objects.erase(itr);
        } else if (neighbors.ny != -1 && (*itr)->y() < bounds.ny) {
            // Copy the object's id and state into transfer buffers
            store_object_state(&idBufferNY, &stateBufferNY, (*itr));
            // Delete the object reference from the local list
            itr = objects.erase(itr);
        } else if (neighbors.py != -1 && (*itr)->y() > bounds.py) {
            // Copy the object's id and state into transfer buffers
            store_object_state(&idBufferPY, &stateBufferPY, (*itr));
            // Delete the object reference from the local list
            itr = objects.erase(itr);
        } else if (neighbors.nz != -1 && (*itr)->z() < bounds.nz) {
            // Copy the object's id and state into transfer buffers
            store_object_state(&idBufferNZ, &stateBufferNZ, (*itr));
            // Delete the object reference from the local list
            itr = objects.erase(itr);
        } else if (neighbors.pz != -1 && (*itr)->z() > bounds.pz) {
            // Copy the object's id and state into transfer buffers
            store_object_state(&idBufferPZ, &stateBufferPZ, (*itr));
            // Delete the object reference from the local list
            itr = objects.erase(itr);
        } else {
            // Advance the iterator
            ++itr;
        }
    }

    // Send object lists
    if (neighbors.nx != -1) send_objects(&idBufferNX, &stateBufferNX, neighbors.nx, nxSendRequestList);
    if (neighbors.px != -1) send_objects(&idBufferPX, &stateBufferPX, neighbors.px, pxSendRequestList);
    if (neighbors.ny != -1) send_objects(&idBufferNY, &stateBufferNY, neighbors.ny, nySendRequestList);
    if (neighbors.py != -1) send_objects(&idBufferPY, &stateBufferPY, neighbors.py, pySendRequestList);
    if (neighbors.nz != -1) send_objects(&idBufferNZ, &stateBufferNZ, neighbors.nz, nzSendRequestList);
    if (neighbors.pz != -1) send_objects(&idBufferPZ, &stateBufferPZ, neighbors.pz, pzSendRequestList);

    // Receive from negative-x neighbors
    if (neighbors.nx != -1) {
        receive_objects(rank-1, &objects);
    }
     // Receive from positive-x neighbors
    if (neighbors.px != -1) {
        receive_objects(rank+1, &objects);
    }
    // Receive from negative-y neighbors
    if (neighbors.ny != -1) {
        receive_objects(neighbors.ny, &objects);
    }
    // Receive from positive-y neighbors
    if (neighbors.py != -1) {
        receive_objects(neighbors.py, &objects);
    }
    // Receive from negative-z neighbors
    if (neighbors.nz != -1) {
        receive_objects(neighbors.nz, &objects);
    }
    // Receive from positive-z neighbors
    if (neighbors.pz != -1) {
        receive_objects(neighbors.pz, &objects);
    }
 
    // Wait for communication to be complete
    if (neighbors.nx != -1) { MPI_Waitall(3, nxSendRequestList, MPI_STATUSES_IGNORE); }
    if (neighbors.px != -1) { MPI_Waitall(3, pxSendRequestList, MPI_STATUSES_IGNORE); }
    if (neighbors.ny != -1) { MPI_Waitall(3, nySendRequestList, MPI_STATUSES_IGNORE); }
    if (neighbors.py != -1) { MPI_Waitall(3, pySendRequestList, MPI_STATUSES_IGNORE); }
    if (neighbors.nz != -1) { MPI_Waitall(3, nzSendRequestList, MPI_STATUSES_IGNORE); }
    if (neighbors.pz != -1) { MPI_Waitall(3, pzSendRequestList, MPI_STATUSES_IGNORE); }
}

void Slice::handle_collision(Object * obj1, Object * obj2) {

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
    for (std::vector<Object *>::iterator itr = objects.begin(); itr != objects.end(); ++itr) {
        for (std::vector<Object *>::iterator ghostItr = ghostObjects.begin(); ghostItr != ghostObjects.end(); ++ghostItr) {
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
    for (std::vector<Object *>::iterator itr = objects.begin(); itr != objects.end(); ++itr) {
        for (std::vector<Object *>::iterator otherItr = itr; otherItr != objects.end(); ++otherItr) {
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
    for (std::vector<Object *>::iterator itr = objects.begin(); itr != objects.end(); ++itr) {
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
        if (neighbors.py == -1 && (*itr)->y() + (*itr)->getRadius() > bounds.py) {
            dp.y = bounds.py - ((*itr)->y() + (*itr)->getRadius());
            dv.y = -1;
        } 
        else if (neighbors.ny == -1 && (*itr)->y() - (*itr)->getRadius() < bounds.ny) {
            dp.y = bounds.ny - ((*itr)->y() - (*itr)->getRadius()); 
            dv.y = -1;
        }
        // z boundaries
        if (neighbors.pz == -1 && (*itr)->z() + (*itr)->getRadius() > bounds.pz) {
            dp.z = bounds.pz - ((*itr)->z() + (*itr)->getRadius());
            dv.z = -1;
        } 
        else if (neighbors.nz == -1 && (*itr)->z() - (*itr)->getRadius() < bounds.nz) {
            dp.z = bounds.nz - ((*itr)->z() - (*itr)->getRadius()); 
            dv.z = -1;
        }
        
        //
        /*
        if (neighbors.pz == -1 && (*itr)->z() + (*itr)->getRadius() > bounds.pz) {
            dp.z = bounds.pz - ((*itr)->z() + (*itr)->getRadius());
            dv.z = -1; 
        } 
        else if (neighbors.nz == -1 && (*itr)->z() - (*itr)->getRadius() < bounds.nz) {
            dp.z = bounds.nz - ((*itr)->z() - (*itr)->getRadius()); 
            dv.z = -1;
        }
        */
        // Update object position
        (*itr)->delta_position(dp);
        (*itr)->multiply_velocity(dv);
    }
}
