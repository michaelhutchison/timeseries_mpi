#include "Mouse.h"

Mouse::Mouse(View * v) {
    mode = 0;
    x = y = 0;
    clickX = clickY = 0;
    view = v;
}
Mouse::~Mouse() {

}

void Mouse::motion(int newx, int newy) {

    x = newx;
    y = newy;

	// Changing horizontal and vertical angle of view
   	if (mode==1)  
   	{
        double newAzimuth = x-clickX;
		double newElevation = y-clickY;
        view->setAngles(newAzimuth, newElevation);
		glutPostRedisplay();
   	}
  	else if (mode==2) 
   	{
        double newAzimuth = x-clickX;
		view->setAzimuth(newAzimuth);
		double newDistance = .01*(y-clickY);
		view->setDistanceFromOrigin(newDistance);
   		view->project();
		glutPostRedisplay();
   	}

}
void Mouse::button(int key,int status,int newx,int newy) {
    //  On button down, set 'mode' and remember location

    x = newx;
    y = newy;

    if (status==GLUT_DOWN)
    {
        if (key==GLUT_LEFT_BUTTON)
        {
            clickX = x - view->getAzimuth();
            clickY = y - view->getElevation();
            mode=1;
        }
        else if (key==GLUT_RIGHT_BUTTON)
        {
            clickX = x - view->getAzimuth();
            clickY = y - (100*view->getDistanceFromOrigin());
            mode=2;
        }
   }
   // set mouse mode to 0 when no mouse button is pressed
   else if (status==GLUT_UP)
      mode = 0;

}
