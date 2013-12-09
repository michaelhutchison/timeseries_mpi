#ifndef RECORDER_H
#define RECORDER_H

#include "vec3.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <stdlib.h>
#include <time.h>

#define RANDOMDOUBLE(min,max) (  ((double)rand() / RAND_MAX) * (max-min) + min  ) 
#ifndef PI
#define PI 3.14159265
#endif
#define UNSIGNEDSHORT sizeof(unsigned short)
#define UNSIGNED sizeof(unsigned)
#define UNSIGNEDLONG sizeof(unsigned long)
#define SHORT sizeof(short)
#define INT sizeof(int)
#define LONG sizeof(long)
#define DOUBLE sizeof(double)

class Object {
public:
    Object();
    ~Object();
    /* getters */
    double x() {return position.x;}
    double y() {return position.y;}
    double z() {return position.z;}
    double a_r() {return rotationAngle;}
    double x_r() {return rotationVector.x;}
    double y_r() {return rotationVector.y;}
    double z_r() {return rotationVector.z;}
    unsigned getID() {return id;}
    Vec3 getPosition() {return position;}
    Vec3 getRotationVector() {return rotationVector;}
    double getRotationAngle() {return rotationAngle;}
    /* setters */
    void setID(unsigned i) {id = i;}
    void setPosition(double x, double y, double z) {position.set(x,y,z);} 
    void setVelocity(double x, double y, double z) {velocity.set(x,y,z);} 
    void setRotationVector(double x, double y, double z) {rotationVector.set(x,y,z);} 
    void setRotationAngle(double a) {rotationAngle = a;}
    void setRotationVelocity(double v) {rotationalVelocity = v;}
    /* Model-specific methods */
    void setAngle(double a) {angle = a;}
    void setRadius(double r) {radius = r;}
    double getAngle() {return angle;}
    double getRadius() {return radius;}
    void advanceFrame();

private:
    unsigned id;
    int objectType;
    /* State data */
    Vec3 position;
    Vec3 rotationVector;
    double rotationAngle;
    /* Transitional data */
    Vec3 velocity;
    double rotationalVelocity;
    /* Model-specific data */
    double radius; // distance from origin
    double angle;  // angle about the z-plane

};


// Class definitions
class ObjectRecord {
public:
    ObjectRecord();
    ObjectRecord(Object obj);
    ~ObjectRecord();
    static const unsigned short objectLengthInBytes = UNSIGNED 
                                             + DOUBLE + DOUBLE + DOUBLE
                                             + DOUBLE + DOUBLE + DOUBLE + DOUBLE;
    Object * getObjectPtr() {return &object;}
    /* IO */
    void write(std::ofstream& fout);
    void read(std::ifstream& fin);
private:
    /* STATE DATA */
    Object object;
};

class FrameState {
public:
    FrameState(unsigned id);
    ~FrameState();
    void set_frame_id(unsigned id) {frameID = id;}
    unsigned get_frame_id()        {return frameID;}
    unsigned get_num_objects()     {return objectRecords.size();}
    Object * get_object(unsigned i);
    void add_object(Object * object);
    void write(std::ofstream& fout);
    void read(std::ifstream& fin);
private:
    std::vector<ObjectRecord *> objectRecords;
    unsigned frameID;
    void write_header(std::ofstream& fout);
    void read_header(std::ifstream& fin);
};


class TimeSeries {
public:
    TimeSeries();
    ~TimeSeries();
    /* Recording */
    void record_frame(Object * objects, unsigned n);
    void write(const char * filename);
    /* Reading */
    void read(const char * filename);
    /* Playback */
    FrameState * get_frame(unsigned i);
    FrameState * get_next_frame();
private:
    std::vector<FrameState *> frames;
    unsigned nFrames;
    unsigned nextFrameIndex;
};

class Model {
public:
    Model(unsigned n);
    ~Model();
    void set_filename(char * newFilename);
    unsigned get_num_objects() {return nObjects;}
    void init();
    void record_frame();
    void finalize_record()  {finalize_record(recordFilename);}
    void finalize_record(const char * filename);
    void advance_frame();
private:
    /* For real time animation */
    unsigned next_object_index;

    char * recordFilename;
    TimeSeries timeseries;
    unsigned nObjects;

    Object * objects;

};


#endif
