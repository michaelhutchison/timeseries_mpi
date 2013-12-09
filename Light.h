#ifndef LIGHT_H
#define LIGHT_H

#include "graphicslib.h"

class Light {
  public:
    Light();
    ~Light();

    void setAmbient(double a) {ambient = a;}
    void setSpecular(double a) {specular = a;}
    void setDiffuse(double a) {diffuse = a;}
    void setColor(double a, double b, double c, double d);

    void draw();
    void setLighting();
  private:
    double x;
    double y;
    double z;
    double ambient;
    double specular;
    double diffuse;
    float lightColor[];
};

#endif
