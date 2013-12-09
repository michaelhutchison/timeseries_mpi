#include "Object_mpi.h"

/******************************
 * Object_mpi
 *  Represents an object's state 
 *  and handles reading and writing the state to/from the recording.
 ******************************/
Object_mpi::Object_mpi() {
    radius = .2;
    rotationAngle = 0;
}
Object_mpi::~Object_mpi() {

}
void Object_mpi::advance_full_step() {
    /*
     * This method defines how the objects in the space change over time.
     * This implementation uses simple explicit Euler update.
     */
    // Update position
    previousPosition = position;
    position += velocity;

    // Update velocity
    velocity += acceleration;

    // Update acceleration
    //  -- not implemented

}
void Object_mpi::advance_half_step() {
    /*
     * This method defines how the objects in the space change over time.
     * This implementation uses simple expclict Euler update.
     */
    // Update position
    previousPosition = position;
    position += .5*velocity;

    // Update velocity
    velocity += .5*acceleration;

    // Update acceleration
    //  -- not implemented

}

    /*****************************************************
     *  OBJECT STATE RECORD
     *  unsigned            object ID
     *  double              x position
     *  double              y position
     *  double              z position
     *  double              x of rotation vector
     *  double              y of rotation vector
     *  double              z of rotation vector
     *  double              angle of rotation
     *****************************************************/
void Object_mpi::get_bytes(char * buf) {
    memcpy(&buf[0], &id, 4);
    memcpy(&buf[4],  &position.x, 8);
    memcpy(&buf[12], &position.y, 8);
    memcpy(&buf[20], &position.z, 8);
    memcpy(&buf[28], &rotationVector.x, 8);
    memcpy(&buf[36], &rotationVector.y, 8);
    memcpy(&buf[44], &rotationVector.z, 8);
    memcpy(&buf[52], &rotationAngle, 8);

}

