#include "World.h"

/* Constructor */
World::World(std::ofstream * fh, double worldSize[3]) {
    fileHandle = fh;

    // initialize world boundaries
    // -- This might be better to initialize outside this class.
    bounds.nx = -worldSize[0] * 0.5;
    bounds.px =  worldSize[0] * 0.5;
    bounds.ny = -worldSize[1] * 0.5;
    bounds.py =  worldSize[1] * 0.5;
    bounds.nz = -worldSize[2] * 0.5;
    bounds.pz =  worldSize[2] * 0.5;

    // Calculate id range
    minID = 0;
    maxID = 4294967290; // 32 bits worth, give or take a few.
    nextObjectID = minID;
    nextFrameID = 0;
}
World::~World() {
    // Release dynamically allocated memory
    for (unsigned i=0; i<objects.size(); i++) 
        delete objects[i];
}
void World::write_frame_header() {
    unsigned short objectLength = Object_mpi::objectLengthInBytes;
    unsigned long frameSizeInBytes = objects.size() * objectLength;

    // write frame size in bytes
    (*fileHandle).write((char*)&frameSizeInBytes, TYPEUNSIGNEDLONG);
    // write frame ID
    (*fileHandle).write((char*)&nextFrameID, TYPEUNSIGNED);
    // Write number of objects
    unsigned nObjects = objects.size();
    (*fileHandle).write((char*)&nObjects, TYPEUNSIGNED);
    // Write length of an object
    (*fileHandle).write((char*)&objectLength, TYPEUNSIGNEDSHORT);


    nextFrameID++;

}
void World::record_frame() {
    // Save the current frame to the animation record

    write_frame_header();
    
    // Create buffer containing all object states
    unsigned totalBytes = objects.size() * Object_mpi::objectLengthInBytes;
    char * objectBuffer = new char[totalBytes];
    memset(objectBuffer, 0, sizeof(objectBuffer));
    for (unsigned i=0; i<objects.size(); i++) {
        objects[i]->get_bytes(&objectBuffer[i*60]);
    }


    (*fileHandle).write(objectBuffer, totalBytes);

}


void World::createObject() {
    // Do nothing if ID range has been used up
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

void World::advance_full_step() {
    /* Advances position and velocity a half-step for each object */
    for (unsigned i=0; i<objects.size(); i++) {
        objects[i]->advance_full_step();
    }   
}
void World::advance_half_step() {
    /* Advances position and velocity a half-step for each object */
    for (unsigned i=0; i<objects.size(); i++) {
        objects[i]->advance_half_step();
    }   
}

void World::handle_collision(Object_mpi * obj1, Object_mpi * obj2) {

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
void World::detect_collisions() {
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
void World::detect_collisions_world_boundaries() {
    
    // NOTE: only handles collisions with world boundaries.
    for (std::vector<Object_mpi *>::iterator itr = objects.begin(); itr != objects.end(); ++itr) {
        Vec3 dp;
        Vec3 dv(1,1,1);

        // x boundaries
        if ((*itr)->x() + (*itr)->getRadius() > bounds.px) {
            dp.x = bounds.px - ((*itr)->x() + (*itr)->getRadius());
            dv.x = -1;
        } 
        else if ((*itr)->x() - (*itr)->getRadius() < bounds.nx) {
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
