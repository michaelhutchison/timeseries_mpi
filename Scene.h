#ifndef SCENE_H
#define SCENE_H

#include "graphicslib.h"
#include "Texture.h"
#include "Light.h"
#include "Model.h"

class Scene {
  public:
    Scene();
    ~Scene();

    void axesOn()       {drawAxesOn = true;}
    void axesOff()      {drawAxesOn = false;}
    void toggleAxes()   {drawAxesOn = !drawAxesOn;}
    void display();
    void initialize();
    void advanceFrame();
private:
    int W;
    int H;
    Texture textures;
    Light light;
    TimeSeries reader;
    Model model;
 
    void drawScene();
    void drawAxes();
    bool drawAxesOn;
    float backgroundColor[4];
    void loadTextures();
    void Cube (double x,double y,double z,
               double xLength, double yLength, double zLength,
               double th,double rx,double ry,double rz);
};

#endif
