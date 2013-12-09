/*
 *  main.cpp performs the following tasks:
 *  1) Sets up GLUT
 *  2) Contains the Scene, View, and Mouse objects 
 *  3) Signals Scene object to update frames at the correct time
 *  4) Detects mouse and keyboard input and delegates to the 
 *     appropriate object for the response.
 */
#include "graphicslib.h"
#include "View.h"
#include "Scene.h"
#include "Mouse.h"

using namespace std;

// Control animation
#define FRAMESPERSECOND 30
#define MSPERFRAME (1000/FRAMESPERSECOND)
int prevTime=0; 	//  Time of last frame update

// Declare objects
Scene scene;
View view;
Mouse mouse(&view);

/*
 *  GLUT calls this routine to display the scene
 */
void display() {
   	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
   	glEnable(GL_DEPTH_TEST);
   	glLoadIdentity();
    // Set up View 
	view.setView();
	// Draw scene
    scene.display();
   	// Write debugging parameters -- (0,0) is lower left corner.
	/*
    glColor3f(0,0,0);
	glWindowPos2i(5,10);
	Print("x: %d  y:%d  mode:%d " , mouse.getx(), mouse.gety(), mouse.getMode());
	*/
    // Render the scene and make it visible
  	glFlush();   	
   	glutSwapBuffers();
}

/*
 * Idle function - controls animation
 */
void idle() {
    int t = glutGet(GLUT_ELAPSED_TIME);
   	if (t-prevTime >= MSPERFRAME) 
   	{
		// set time of last frame change
		prevTime=t;
        // Update animated elements
        scene.advanceFrame(); 
		// Redraw the scene
		glutPostRedisplay();
   	}
}

/*
 *  GLUT calls this routine when an arrow key, function key, or PageUp/PageDown is pressed
 */
void special(int key,int x,int y) {
 	//  Arrow keys 
   	if (key == GLUT_KEY_RIGHT)  { }
   	else if (key == GLUT_KEY_LEFT)   { }
   	else if (key == GLUT_KEY_UP)     { }
   	else if (key == GLUT_KEY_DOWN)   { }
    //  PageUp and PageDown keys 
    else if (key == GLUT_KEY_PAGE_UP)
      view.setDimension(view.getDimension() - 0.1);
    else if (key == GLUT_KEY_PAGE_DOWN )
      view.setDimension(view.getDimension() + 0.1);
    // Other special keys
    else if (key == GLUT_KEY_HOME)    { }
    else if (key == GLUT_KEY_END)     { }
    else if (key == GLUT_KEY_INSERT)  { }
    // F-keys 
   	else if (key == GLUT_KEY_F1)    { }
   	else if (key == GLUT_KEY_F2)    { }
   	else if (key == GLUT_KEY_F3)    { }
   	else if (key == GLUT_KEY_F4)    { }
   	else if (key == GLUT_KEY_F5)    { }
   	else if (key == GLUT_KEY_F6)    { }
   	else if (key == GLUT_KEY_F7)    { }
   	else if (key == GLUT_KEY_F8)    { }
   	else if (key == GLUT_KEY_F9)    { }
   	else if (key == GLUT_KEY_F10)   { }
   	else if (key == GLUT_KEY_F11)   { }
   	else if (key == GLUT_KEY_F12)   { }
   	//  Update projection and redraw
   	view.project();
   	glutPostRedisplay();
}

/*
 *  GLUT calls this routine when a key is pressed
 *  NOTE: escape, backspace, and delete keys are generated as ASCII characters.
 */
void key(unsigned char ch,int x,int y) {
	// escape to quit program
   	if (ch == 27) 	
		exit(0);
   	// toggle axes
   	else if (ch =='x' || ch=='X')
        scene.toggleAxes();
   	// Update projection and redraw
   	view.project();
   	glutPostRedisplay();
}

/*
 *  GLUT calls this routine when a mouse is moved while a button is pressed
 */
void mouseButtonMotion(int x,int y) {
    mouse.motion(x,y);
}
/*
 *  GLUT calls this routine when a mouse button is pressed or released
 */
void mouseButton(int key,int status,int x,int y) {
    mouse.button(key, status, x,  y);
}
/*
 *  GLUT calls this routine when a mouse is moved while no button is pressed
 */
void mouseMotion(int x,int y) {
    mouse.motion(x,y);
}
/*
 * GLUT calls this routine when a window is reshaped. 
 */
void reshape(int width, int height) {
    view.reshape(width,height);
}
/*
 *  Start up and configure GLUT
 */
int main(int argc,char* argv[]) {
   	// Initialize GLUT
   	glutInit(&argc,argv);
	// Request double buffered, true color window with Z buffering at 600x600
	glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
	glutInitWindowSize(800,800);
	glutCreateWindow("Timeseries Playback");
	// Set glut callbacks
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutSpecialFunc(special);
	glutKeyboardFunc(key);
	glutIdleFunc(idle);
	glutMouseFunc(mouseButton);
	glutMotionFunc(mouseButtonMotion);
    glutPassiveMotionFunc(mouseMotion);
    // Initialize the scene
    scene.initialize();
    // Initialize the View
    view.setAzimuth(102);       //  Azimuth of view angle
    view.setElevation(45);      //  Elevation of view angle
    view.setFieldOfView(35);    //  Field of view (for perspective)
    view.setAspectRatio(1);     //  Aspect ratio
    view.setDimension(14.0);    //  Size of world
	// Pass control to GLUT so it can interact with the user
	glutMainLoop();
	return 0;
}
