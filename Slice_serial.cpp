#include "Slice_serial.h"

/* Constructor */
Slice_serial::Slice_serial(unsigned r, unsigned s, std::ofstream * fh, double worldSize[3], unsigned sliceRows[3]) {
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

    // Initialize neighbors pointers
    // -- value of  indicates no neighbor at that boundary
    neighbors.nx = 0;
    neighbors.px = 0;
    neighbors.ny = 0;
    neighbors.py = 0;
    neighbors.nz = 0;
    neighbors.pz = 0;
    neighbors_rank.nx = -1;
    neighbors_rank.px = -1;
    neighbors_rank.ny = -1;
    neighbors_rank.py = -1;
    neighbors_rank.nz = -1;
    neighbors_rank.pz = -1;
    if (rank > minInRowX) neighbors_rank.nx = rank-1;
    if (rank < maxInRowX) neighbors_rank.px = rank+1;
    if (rank > minInRowY) neighbors_rank.ny = rank-slicesX;
    if (rank < maxInRowY) neighbors_rank.py = rank+slicesX;
    if (rank > minInRowZ) neighbors_rank.nz = rank-(slicesX*slicesY);
    if (rank < maxInRowZ) neighbors_rank.pz = rank+(slicesX*slicesY);
    
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

Slice_serial::~Slice_serial() {
    // Release dynamically allocated memory
    for (unsigned i=0; i<objects.size(); i++) 
        delete objects[i];
}
void Slice_serial::set_neighbors(std::vector<Slice_serial *> slices) {
    // In the constructor we can determine the rank of any neighbors.
    // We cannot get a pointer to neighboring Slice objects until
    // all Slices have been created. 
    // After this action, we probably no longer need to know their ranks.
    if (neighbors_rank.nx != -1) neighbors.nx = slices[neighbors_rank.nx]; 
    if (neighbors_rank.px != -1) neighbors.px = slices[neighbors_rank.px]; 
    if (neighbors_rank.ny != -1) neighbors.ny = slices[neighbors_rank.ny]; 
    if (neighbors_rank.py != -1) neighbors.py = slices[neighbors_rank.py]; 
    if (neighbors_rank.nz != -1) neighbors.nz = slices[neighbors_rank.nz]; 
    if (neighbors_rank.pz != -1) neighbors.pz = slices[neighbors_rank.pz]; 
}
void Slice_serial::write_frame_header() {
    unsigned short objectLength = Object::objectLengthInBytes;
    unsigned long frameSizeInBytes = nTotalObjects * objectLength;

    // write frame size in bytes
    (*fileHandle).write((char*)&frameSizeInBytes, TYPEUNSIGNEDLONG);
    // write frame ID
    (*fileHandle).write((char*)&nextFrameID, TYPEUNSIGNED);
    // Write number of objects
    unsigned nObjects = objects.size();
    (*fileHandle).write((char*)&nTotalObjects, TYPEUNSIGNED);
    // Write length of an object
    (*fileHandle).write((char*)&objectLength, TYPEUNSIGNEDSHORT);

    nextFrameID++;
}
void Slice_serial::record_frame() {
    // Save the current frame to the animation record
    if (rank == 0) write_frame_header();
    
    // Create buffer containing all object states
    unsigned totalBytes = objects.size() * Object::objectLengthInBytes;
    char * objectBuffer = new char[totalBytes];
    memset(objectBuffer, 0, sizeof(objectBuffer));
    for (unsigned i=0; i<objects.size(); i++) {
        objects[i]->get_bytes(&objectBuffer[i*60]);
    }


    (*fileHandle).write(objectBuffer, totalBytes);

}


void Slice_serial::createObject() {
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

void Slice_serial::advance_full_step() {
    /* Advances position and velocity a half-step for each object */
    for (unsigned i=0; i<objects.size(); i++) {
        objects[i]->advance_full_step();
    }   
}
void Slice_serial::advance_half_step() {
    /* Advances position and velocity a half-step for each object */
    for (unsigned i=0; i<objects.size(); i++) {
        objects[i]->advance_half_step();
    }   
}
void Slice_serial::clear_ghost_objects() {
    ghostObjects.clear();
}
void Slice_serial::exchange_ghost_objects() {
    // Any objects which have crossed the inner boundary with another slice
    // are transferred to that slice as ghost objects.
   
    // Examine each object
    for (std::vector<Object *>::iterator itr = objects.begin(); itr != objects.end(); ++itr) {
        if (neighbors.nx != 0 && (*itr)->x() < innerBounds.nx) {
            neighbors.nx->addGhostObject(*itr);
        } 
        if (neighbors.px != 0 && (*itr)->x() > innerBounds.px) {
            neighbors.px->addGhostObject(*itr);
        }
        if (neighbors.ny != 0 && (*itr)->y() < innerBounds.ny) {
            neighbors.ny->addGhostObject(*itr);
        }
        if (neighbors.py != 0 && (*itr)->y() > innerBounds.py) {
            neighbors.py->addGhostObject(*itr);
        } 
        if (neighbors.nz != 0 && (*itr)->z() < innerBounds.nz) {
            neighbors.nz->addGhostObject(*itr);
        }
        if (neighbors.pz != 0 && (*itr)->z() > innerBounds.pz) {
            neighbors.pz->addGhostObject(*itr);
        } 
 
    }
}

void Slice_serial::exchange_objects() {
    // Any objects which have crossed a boundary with another slice
    // are transferred to that slice. 

    // Examine each object
    for (std::vector<Object *>::iterator itr = objects.begin(); itr != objects.end(); ) {
        if (neighbors.nx != 0 && (*itr)->x() < bounds.nx) {
            neighbors.nx->addObject(*itr);
            // Delete the object reference from the local list
            itr = objects.erase(itr);
        } else if (neighbors.px != 0 && (*itr)->x() > bounds.px) {
            neighbors.px->addObject(*itr);
            // Delete the object reference from the local list
            itr = objects.erase(itr);
        } else if (neighbors.ny != 0 && (*itr)->y() < bounds.ny) {
            neighbors.ny->addObject(*itr);
            // Delete the object reference from the local list
            itr = objects.erase(itr);
        } else if (neighbors.py != 0 && (*itr)->y() > bounds.py) {
            neighbors.py->addObject(*itr);
            // Delete the object reference from the local list
            itr = objects.erase(itr);
        } else if (neighbors.nz != 0 && (*itr)->z() < bounds.nz) {
            neighbors.nz->addObject(*itr);
            // Delete the object reference from the local list
            itr = objects.erase(itr);
        } else if (neighbors.pz != 0 && (*itr)->z() > bounds.pz) {
            neighbors.pz->addObject(*itr);
            // Delete the object reference from the local list
            itr = objects.erase(itr);
        } else {
            // Advance the iterator
            ++itr;
        }
    }
}

void Slice_serial::handle_collision(Object * obj1, Object * obj2) {

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
void Slice_serial::detect_collisions() {
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
void Slice_serial::detect_collisions_world_boundaries() {
    
    // NOTE: only handles collisions with world boundaries.
    for (std::vector<Object *>::iterator itr = objects.begin(); itr != objects.end(); ++itr) {
        Vec3 dp;
        Vec3 dv(1,1,1);

        // x boundaries
        if (neighbors.px == 0 && (*itr)->x() + (*itr)->getRadius() > bounds.px) {
            dp.x = bounds.px - ((*itr)->x() + (*itr)->getRadius());
            dv.x = -1;
        } 
        else if (neighbors.nx == 0 && (*itr)->x() - (*itr)->getRadius() < bounds.nx) {
            dp.x = bounds.nx - ((*itr)->x() - (*itr)->getRadius()); 
            dv.x = -1;
        }
        // y boundaries
        if (neighbors.py == 0 && (*itr)->y() + (*itr)->getRadius() > bounds.py) {
            dp.y = bounds.py - ((*itr)->y() + (*itr)->getRadius());
            dv.y = -1;
        } 
        else if (neighbors.ny == 0 && (*itr)->y() - (*itr)->getRadius() < bounds.ny) {
            dp.y = bounds.ny - ((*itr)->y() - (*itr)->getRadius()); 
            dv.y = -1;
        }
        // z boundaries
        if (neighbors.pz == 0 && (*itr)->z() + (*itr)->getRadius() > bounds.pz) {
            dp.z = bounds.pz - ((*itr)->z() + (*itr)->getRadius());
            dv.z = -1;
        } 
        else if (neighbors.nz == 0 && (*itr)->z() - (*itr)->getRadius() < bounds.nz) {
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
