#include "View.h"

View::View() {
    azimuthAngle=102;   //  Azimuth of view angle
    elevationAngle=45;  //  Elevation of view angle
    fieldOfView=35;     //  Field of view (for perspective)
    aspectRatio=1;      //  Aspect ratio
    dimension=14.0;     //  Size of world
    distanceFromOrigin = 1.8;

    setEyePosition();
    vx = 0;
    vy = 0;
    vz = 0;

}
View::~View() {

}
/*
 *  Configure view mode
 */
void View::setView() {
    gluLookAt(  ex,ey,ez, 
                vx,vy,vz, 
                0,Cos(elevationAngle),0);
}
void View::project() {
   //  Tell OpenGL we want to manipulate the projection matrix
   glMatrixMode(GL_PROJECTION);
   //  Undo previous transformations
   glLoadIdentity();
   //  Perspective transformation
   if (fieldOfView)
      gluPerspective(fieldOfView,aspectRatio,dimension/4,4*dimension);
   //  Orthogonal transformation
   else
      glOrtho(-aspectRatio*dimension,
               aspectRatio*dimension,
              -dimension,
              +dimension,
              -dimension,
              +dimension);
   //  Switch to manipulating the model matrix
   glMatrixMode(GL_MODELVIEW);
   //  Undo previous transformations
   glLoadIdentity();

}
/*
 *  GLUT calls this routine when the window is resized
 */
void View::reshape(int width,int height) {
   	//  Ratio of the width to the height of the window
   	aspectRatio = (height>0) ? (double)width/height : 1;
   	//  Set the viewport to the entire window
   	glViewport(0,0, width,height);
   	//  Set projection
    project();
}
void View::setEyePosition() {
    ex = -distanceFromOrigin*dimension*Sin(azimuthAngle)*Cos(elevationAngle);
    ey = +distanceFromOrigin*dimension                  *Sin(elevationAngle);
    ez = +distanceFromOrigin*dimension*Cos(azimuthAngle)*Cos(elevationAngle);
 
}
void View::setDistanceFromOrigin(float d) {
    if (d > 0.5){ 
        distanceFromOrigin = d;
        fieldOfView = (10.94*distanceFromOrigin)+15.624;
        setEyePosition();
    }
}
void View::setAzimuth(int az) {
    azimuthAngle = az%360;
    setEyePosition();
}
void View::setElevation(int el) {
    elevationAngle = el%360;
    setEyePosition();
}
void View::setAngles(int az, int el){
    azimuthAngle = az%360;
    elevationAngle = el%360;
    setEyePosition();
}

