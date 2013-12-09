#include "Model.h"

/******************************
 * ObjectRecord
 *  Writes the state of a single object in the scene 
 *  to the timeseries record.
 ******************************/
ObjectRecord::ObjectRecord() {

}
ObjectRecord::ObjectRecord(Object obj) {
    object = obj;
}
ObjectRecord::~ObjectRecord() {

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
void ObjectRecord::write(std::ofstream& fout) {
    // Copy data from object
    unsigned id = object.getID();
    Vec3 position = object.getPosition();
    Vec3 rotation = object.getRotationVector();
    double angle = object.getRotationAngle();
    // Write object ID
    fout.write((char*)&id, UNSIGNED);
    // Write position
    fout.write((char*)&position.x,DOUBLE);
    fout.write((char*)&position.y,DOUBLE);
    fout.write((char*)&position.z,DOUBLE);
    // Write orientation 
    fout.write((char*)&rotation.x,DOUBLE);
    fout.write((char*)&rotation.y,DOUBLE);
    fout.write((char*)&rotation.z,DOUBLE);
    fout.write((char*)&angle,DOUBLE);
    /* 
    std::cout << "   object ID " << objectID << std::endl;
    std::cout << "   x pos " << xPos << std::endl;
    std::cout << "   y pos " << yPos << std::endl;
    std::cout << "   z pos " << zPos << std::endl;
    std::cout << "   x rot " << xRot << std::endl;
    std::cout << "   y rot " << yRot << std::endl;
    std::cout << "   z rot " << zRot << std::endl;
    std::cout << "   a rot " << aRot << std::endl;
    */
}
void ObjectRecord::read(std::ifstream& fin) {
    // Create local buffers for incoming data
    unsigned id;
    Vec3 position;
    Vec3 rotation;
    double angle;
    // Read object ID
    fin.read((char*)&id, UNSIGNED);
    //std::cout << "Object id: " << id << std::endl;
    // Read position
    fin.read((char*)&position.x, DOUBLE);
    fin.read((char*)&position.y, DOUBLE);
    fin.read((char*)&position.z, DOUBLE);
    //std::cout << "Object position: " << position.x << "  " << position.y << "  " << position.z << std::endl;
    // Read orientation
    fin.read((char*)&rotation.x, DOUBLE);
    fin.read((char*)&rotation.y, DOUBLE);
    fin.read((char*)&rotation.z, DOUBLE);
    fin.read((char*)&angle, DOUBLE);
    //std::cout << "orientation: " << rotation.x << "  " << rotation.y << "  " << rotation.z << std::endl;
    //std::cout << "angle: " << angle << std::endl;

    // Copy local data to object state
    object.setID(id);
    object.setPosition(position.x, position.y, position.z);
    object.setRotationVector(rotation.x, rotation.y, rotation.z);
    object.setRotationAngle(angle);
    /*
    std::cout << "   object ID " << objectID << std::endl;
    std::cout << "   x pos " << xPos << std::endl;
    std::cout << "   y pos " << yPos << std::endl;
    std::cout << "   z pos " << zPos << std::endl;
    std::cout << "   x rot " << xRot << std::endl;
    std::cout << "   y rot " << yRot << std::endl;
    std::cout << "   z rot " << zRot << std::endl;
    std::cout << "   a rot " << aRot << std::endl;
    */
}

/******************************
 * FrameState
 *  Stores position and orientation (ObjectState)
 *  for a collection of objects in a single 
 *  time slice of the scene.
 ******************************/
FrameState::FrameState(unsigned id) {
    frameID = id;
}
FrameState::~FrameState() {
    // Free dynamically allocated memory
    for (unsigned i=0; i < objectRecords.size(); i++) {
        delete objectRecords[i];
    }
}
Object * FrameState::get_object(unsigned i) {
    if (i >= objectRecords.size()) return 0;
    Object * o = objectRecords[i]->getObjectPtr(); 
    return o;
}
void FrameState::add_object(Object * object) {
    ObjectRecord * o = new ObjectRecord(*object);
    objectRecords.push_back(o);
}
void FrameState::write(std::ofstream& fout) {
    //std::cout << "FRAME " << frameID << std::endl;
    write_header(fout);
    for (unsigned i=0; i < objectRecords.size(); i++) {
        objectRecords[i]->write(fout);
    }
}
void FrameState::read(std::ifstream& fin){
    //std::cout << "FRAME IN" << std::endl;
    read_header(fin);
    for (unsigned i=0; i<objectRecords.size(); i++){
        ObjectRecord * o = new ObjectRecord();
        objectRecords[i] = o;
        objectRecords[i]->read(fin);
    }
}
    /*****************************************************
     *  FRAME STATE HEADER
     *  unsigned long       frame size in bytes
     *  unsigned int        frame ID
     *  unsigned int        number of objects in frame
     *  unsigned short      length of each object in bytes
     *****************************************************/
void FrameState::write_header(std::ofstream& fout) {
   // Calculate data sizes
    unsigned short singleObjectLengthInBytes = ObjectRecord::objectLengthInBytes;
    unsigned long allObjectsLengthInBytes = singleObjectLengthInBytes * objectRecords.size();
    unsigned long frameLengthInBytes = UNSIGNEDLONG 
                                      + UNSIGNED 
                                      + UNSIGNED
                                      + UNSIGNEDSHORT
                                      + allObjectsLengthInBytes;
    unsigned nObjects = objectRecords.size();

    
    // Write frame size in bytes
    fout.write((char*)&frameLengthInBytes,UNSIGNEDLONG);
    // Write frame ID
    fout.write((char*)&frameID,UNSIGNED);
    // Write number of objects
    fout.write((char*)&nObjects,UNSIGNED);
    // Write length of an object
    fout.write((char*)&singleObjectLengthInBytes,UNSIGNEDSHORT);
    /*
    std::cout << " frame len " << frameLengthInBytes << std::endl;
    std::cout << " frame ID " << frameID << std::endl;
    std::cout << " nObjects " << nObjects << std::endl;
    std::cout << " object len " << singleObjectLengthInBytes << std::endl;
    */
}
void FrameState::read_header(std::ifstream& fin) {
    // Read frame size in bytes
    unsigned long frameLengthInBytes;
    fin.read((char*)&frameLengthInBytes,UNSIGNEDLONG);
    //std::cout << "Frame length: " << frameLengthInBytes << std::endl;
    // Read frame ID
    fin.read((char*)&frameID,UNSIGNED);
    std::cout << "Frame id: " << frameID << std::endl;
    // Read number of objects
    unsigned nObjects;
    fin.read((char*)&nObjects,UNSIGNED);
    std::cout << "# objects: " << nObjects << std::endl;
    // Read length of an object
    unsigned short singleObjectLengthInBytes;
    fin.read((char*)&singleObjectLengthInBytes,UNSIGNEDSHORT);
    //std::cout << "object length: " << singleObjectLengthInBytes << std::endl;
    /*
    std::cout << " frame len " << frameLengthInBytes << std::endl;
    std::cout << " frame ID " << frameID << std::endl;
    std::cout << " nObjects " << nObjects << std::endl;
    std::cout << " object len " << singleObjectLengthInBytes << std::endl;
    */
    // Clean and re-size objects vector
    for (unsigned i=0; i<objectRecords.size(); i++) {
        delete objectRecords[i];
    }
    objectRecords.resize(nObjects);
}


/******************************
 * TimeSeries
 *  Stores a series of precomputed
 *  frames composing an animation.
 ******************************/
TimeSeries::TimeSeries() {
    nFrames = 0;
    nextFrameIndex = 0;
}
TimeSeries::~TimeSeries() {
    // Release dynamically allocated memory
    for (unsigned i=0; i < frames.size(); i++) {
        delete frames[i];
    }
}
FrameState * TimeSeries::get_frame(unsigned i) {
    if (i >= frames.size()) return 0;
    return frames[i];
}
FrameState * TimeSeries::get_next_frame() {
    FrameState * nextFrame = get_frame(nextFrameIndex);
    nextFrameIndex++;
    if (nextFrameIndex >= frames.size()) nextFrameIndex = 0;
    return nextFrame;
}

void TimeSeries::record_frame(Object * objects, unsigned n) {;
    // Create a new frame
    FrameState * frame = new FrameState(nFrames);
    frames.push_back(frame);
    nFrames++;

    // Record each object in the frame
    for (unsigned i=0; i<n; i++) {
        Object * obj_ptr = &objects[i];
        frame->add_object(obj_ptr);
    }
}
    /*****************************************************
     *  TIME SERIES FILE HEADER
     *  unsigned short      header length in bytes
     *  unsigned long       number of frames in file
     *****************************************************/
void TimeSeries::write(const char * filename) {
    std::cout << "Writing timeseries data" << std::endl;
    // Open file for writing
    std::ofstream fout;
    fout.open(filename, std::ios::out|std::ios::binary);
    // Write header length
    unsigned short headerLengthInBytes = UNSIGNEDSHORT 
                                       + UNSIGNED;
    fout.write((char*)&headerLengthInBytes,UNSIGNEDSHORT);
    // Write number of frames
    fout.write((char*)&nFrames,UNSIGNED);
    /*
    std::cout << " header len " << headerLengthInBytes << std::endl;
    std::cout << " nFrames " << nFrames << std::endl;
    */
    // Write each frame
    for (unsigned i=0; i<frames.size(); i++) {
        frames[i]->write(fout);
    }
    // Cleanup
    fout.close();
}
void TimeSeries::read(const char * filename) {
    std::cout << "Loading timeseries data" << std::endl;
    // Reset: Clear frames and initialize metadata
    for (unsigned i=0; i < frames.size(); i++) {
        delete frames[i];
    }
    nFrames = 0;
    nextFrameIndex = 0;
    // Open file for reading
    std::ifstream fin;
    fin.open(filename, std::ios::in|std::ios::binary);
    // Read header length in bytes
    unsigned short headerLengthInBytes;
    fin.read((char*)&headerLengthInBytes, UNSIGNEDSHORT);
    std::cout << "Header length in bytes: " << headerLengthInBytes << std::endl;
    // Read number of frames
    fin.read((char*)&nFrames, UNSIGNED);
    std::cout << "Number of frames: " << nFrames << std::endl;
    /*
    std::cout << " header len " << headerLengthInBytes << std::endl;
    std::cout << " nFrames " << nFrames << std::endl;
    */
    // Read each frame
    for (unsigned i=0; i<nFrames; i++){
        // Create a new frame
        FrameState * frame = new FrameState(i);
        frames.push_back(frame);
        frame->read(fin);
    }
    // Cleanup
    fin.close();
}


/******************************
 * Object
 *  
 ******************************/
Object::Object() {
    angle = 0;
    radius = 0;
}
Object::~Object() {

}
void Object::advanceFrame() {
    /*
     * This method defines how the objects in the space change over time.
     *   objects rotate around the origin at a speed inversely proportional 
     *      to their distance from the origin.
     *   objects spin on a unique rotational axis at a steady angular velocity.
     */
 
    // Update angle about the origin
    double angularVelocity = 0.05/radius;
    if (angularVelocity <= 0) angularVelocity = 0.01;
    angle += angularVelocity;
    while (angle >= 2*PI) angle -= 2*PI;

    // Update position - fixed distance from origin at the new angle
    position.x = radius * cos(angle);
    position.y = radius * sin(angle);

    // Update rotation angle - rotation vector does not change
    rotationAngle += rotationalVelocity;
    while (rotationAngle >= 2*PI) rotationAngle -= 2*PI;
}

/******************************
 * Model
 *  
 ******************************/
/* Constructor */
Model::Model(unsigned n) {
    nObjects = n;
    recordFilename = new char[101];
    strncpy(recordFilename, "animation.rcd",100);
    next_object_index = 0;
    // Allocate memory
    objects = new Object[nObjects];

}
/* Destructor */
Model::~Model() {
    // Release dynamically allocated memory
    delete[] objects;
    delete[] recordFilename;
}

void Model::set_filename(char * newFilename) {
    strncpy(recordFilename, newFilename, 100);
}
void Model::init() {
    // randomize time
    srand (time(NULL));
    // Generate initial state
    for (unsigned i=0; i<nObjects; i++) {
        // Initialize object i
        double angle = RANDOMDOUBLE(0, 2*PI);
        objects[i].setAngle(angle);
        double radius = RANDOMDOUBLE(0.1, 5.0);
        objects[i].setRadius(radius);
        objects[i].setPosition(radius*cos(angle),
                               radius*sin(angle),
                               RANDOMDOUBLE(-3.0,3.0) );
        objects[i].setRotationVector(RANDOMDOUBLE(-1.0,1.0),
                               RANDOMDOUBLE(-1.0,1.0),
                               RANDOMDOUBLE(-1.0,1.0));
        objects[i].setRotationAngle(RANDOMDOUBLE(0,2*PI));
        objects[i].setRotationVelocity(RANDOMDOUBLE(-5,5));
    }
}
void Model::record_frame() {
    // Save the current frame to the animation record
    timeseries.record_frame(objects, nObjects);

}
void Model::finalize_record(const char * filename) {
    // Write all saved timeseries frames to file
    timeseries.write(filename);
}
void Model::advance_frame() {
    for (unsigned i=0; i<nObjects; i++) {
        objects[i].advanceFrame();
    }
    next_object_index = 0;
}
