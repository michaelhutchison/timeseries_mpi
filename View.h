#ifndef VIEW_H
#define VIEW_H

#include "graphicslib.h"
#include "vec3.h"

class View {

  public:
    View();
    ~View();

    void setAngles(int az, int el);
    void setAzimuth(int az);
    void setElevation(int el);
    void setFieldOfView(int fov) {fieldOfView = fov;}
    void setAspectRatio(int asp) {aspectRatio = asp;}
    void setDimension(int dim) {dimension = dim;}
    void setDistanceFromOrigin(float d); 
    
    int getAzimuth() {return azimuthAngle;}
    int getElevation() {return elevationAngle;}
    int getFieldOfView() {return fieldOfView;}
    int getAspectRatio() {return aspectRatio;}
    int getDimension() {return dimension;}
    int getDistanceFromOrigin() {return distanceFromOrigin;}

    void reshape(int width,int height);
    void setView();
    void project();
  private:
    void setEyePosition();

    int azimuthAngle;
    int elevationAngle;
    int fieldOfView;
    int aspectRatio;
    int dimension;
    double distanceFromOrigin;

    double ex,ey,ez;  // position of viewer
    double vx,vy,vz;  // look-at location

};

#endif
