#ifndef SCENE_H
#define SCENE_H

#include "graphicslib.h"
#include "Texture.h"
#include "Light.h"
#include "Model.h"
#include "vec3.h"

class Scene {
  public:
    Scene();
    ~Scene();

    void axesOn()       {drawAxesOn = true;}
    void axesOff()      {drawAxesOn = false;}
    void toggleAxes()   {drawAxesOn = !drawAxesOn;}
    void display();
    void initialize(char * file);
    void advanceFrame();
    void setWorldSize(double size[3]);
private:
    int W;
    int H;
    double worldX,worldY,worldZ; 
    Texture textures;
    Light light;
    TimeSeries reader;
    Model model;
    char filename[100];
 
    void drawScene();
    void drawAxes() {drawAxes(0,0,0);}
    void drawAxes(double x, double y, double z);
    bool drawAxesOn;
    float backgroundColor[4];
    void loadTextures();

    Vec3 colors[16];
    void WireframeCube (double x,double y,double z,
               double xLength, double yLength, double zLength,
               double th,double rx,double ry,double rz);
 
    void Cube (double x,double y,double z,
               double xLength, double yLength, double zLength,
               double th,double rx,double ry,double rz);
    void Sphere(double x, double y, double z, double radius);
};

#endif
