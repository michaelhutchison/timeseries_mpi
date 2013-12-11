#pragma once

#include "vec3.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <stdlib.h>
#include <time.h>
#include <cmath>

#define TYPEUNSIGNEDSHORT sizeof(unsigned short)
#define TYPEUNSIGNED sizeof(unsigned)
#define TYPEUNSIGNEDLONG sizeof(unsigned long)
#define TYPESHORT sizeof(short)
#define TYPEINT sizeof(int)
#define TYPELONG sizeof(long)
#define TYPEDOUBLE sizeof(double)

class Object_mpi {
public:
    Object_mpi();
    ~Object_mpi();
    /* getters */
    double x() {return position.x;}
    double y() {return position.y;}
    double z() {return position.z;}
    unsigned getID() {return id;}
    Vec3 getPosition() {return position;}
    Vec3 getVelocity() {return velocity;}
    Vec3 getRotationVector() {return rotationVector;}
    double getRotationAngle() {return rotationAngle;}
    double getRadius() {return radius;}

    /* setters */
    void setID(unsigned i) {id = i;}
    void setPosition(double x, double y, double z) {position.set(x,y,z);} 
    void setPosition(Vec3 p) {position = p;}
    void setVelocity(double x, double y, double z) {velocity.set(x,y,z);} 
    void setVelocity(Vec3 v) {velocity = v;}
    void setRotationVector(double x, double y, double z) {rotationVector.set(x,y,z);}
    void setRotationVector(Vec3 v) {rotationVector = v;}
    void setRotationAngle(double a) {rotationAngle = a;}
    void setRotationVelocity(double v) {rotationVelocity = v;}
    /* Transition methods */
    void delta_position(Vec3 dp) {position+=dp;}
    void multiply_velocity(Vec3 dv) {velocity.x*=dv.x; velocity.y*=dv.y; velocity.z*=dv.z;}
    void advance_full_step();
    void advance_half_step();
    /* IO */
    void get_bytes(char * buf);
    static const unsigned short objectLengthInBytes = TYPEUNSIGNED + TYPEDOUBLE + TYPEDOUBLE + TYPEDOUBLE + TYPEDOUBLE + TYPEDOUBLE + TYPEDOUBLE + TYPEDOUBLE;
 
private:
    unsigned id;
    double radius;
    /* State data */
    Vec3 position;
    Vec3 previousPosition;
    Vec3 rotationVector;
    double rotationAngle;
    /* Transitional data */
    Vec3 velocity;
    Vec3 acceleration;
    double rotationVelocity;
};




