#include "Light.h"


Light::Light() {
    ambient=0.1;
    specular=.2;
    diffuse=.4;
    x = 4;
    y = 3;
    z = 0;
    setColor(1.0,1.0,1.0,1.0);

}

Light::~Light() {

}
void Light::draw() {

}
void Light::setColor(double a, double b, double c, double d) {
    lightColor[0] = a;
    lightColor[1] = b;
    lightColor[2] = c;
    lightColor[3] = d;
}
void Light::setLighting() {
    //  Translate intensity to color vectors
    float Ambient[]   = {ambient,ambient,ambient,1.0};
    float Diffuse[]   = {diffuse,diffuse,diffuse,1.0};
    float Specular[]  = {specular,specular,specular,1.0};
	float Position[] = {x,y,z,1.0}; 
   	//  Enable lighting with normalization
    glEnable(GL_NORMALIZE);
    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHTING);
    //  Location of viewer for specular calculations
    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,1);
    //  Enable light 0
    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0,GL_AMBIENT ,Ambient);
    glLightfv(GL_LIGHT0,GL_DIFFUSE ,Diffuse);
    glLightfv(GL_LIGHT0,GL_SPECULAR,Specular);
    glLightfv(GL_LIGHT0,GL_POSITION,Position);


}
