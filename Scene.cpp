#include "Scene.h"


Scene::Scene() : model(1000) {
    axesOn();
    W = 600;
    H = 600;
    backgroundColor[0] = 1.0;
    backgroundColor[1] = 1.0;
    backgroundColor[2] = 1.0;
    backgroundColor[3] = 1.0;
}
Scene::~Scene() {

}
void Scene::initialize() {
    loadTextures();
    reader.read("timeseries.rcd");
    model.init();
}
void Scene::loadTextures() {
    textures.loadTexture("textures/metal2.bmp","metal");
}
void Scene::advanceFrame() {
    //model.advance_frame();
}

void Scene::display() {
    // Set lighting
    light.setLighting();
    // Set background color
    glClearColor(backgroundColor[0],
                 backgroundColor[1],
                 backgroundColor[2],
                 backgroundColor[3]);
 	glColor3f(1,1,1);
    if (drawAxesOn) drawAxes();

    // Draw real time scene
    /*
    Thing thing;
    unsigned long nObjects = model.get_num_objects();
    Vec3 * p = 0;
    Vec3 * r = 0;
    double * a = 0;

    for (unsigned long i=0; i<nObjects; i++) {
        model.get_next_object_state(&p,&r,&a);
        if (p != 0) {
            thing.Cube(p->x, p->y, p->z,
                       0.4,0.4,0.4,
                       *a,p->x,p->y,p->z);
        }
        else { std::cout << "its zero" << std::endl;}
    }
    */

    // Draw recorded scene
    FrameState * frame = reader.get_next_frame();
    if (frame != 0) {
        unsigned long nObjects = frame->get_num_objects();
        for (unsigned long i=0; i<nObjects; i++) {
            Object * o = frame->get_object(i);
            if (o != 0) {
                // Draw an object from the animation record file
                Cube(o->x(), o->y(), o->z(),
                           0.4,0.4,0.4,
                           o->a_r(),o->x_r(),o->y_r(),o->z_r());
            }
        }
    }
}
/*
Draws X Y Z axes
*/
void Scene::drawAxes()
{
	double len=1;
 	glColor3f(1,1,1);
    glBegin(GL_LINES);
    glVertex3d(0.0,0.0,0.0);
    glVertex3d(len,0.0,0.0);
    glVertex3d(0.0,0.0,0.0);
    glVertex3d(0.0,len,0.0);
    glVertex3d(0.0,0.0,0.0);
    glVertex3d(0.0,0.0,len);
    glEnd();
    //  Label axes
    glRasterPos3d(len,0.0,0.0);
    Print("X");
    glRasterPos3d(0.0,len,0.0);
    Print("Y");
    glRasterPos3d(0.0,0.0,len);
    Print("Z");
}
/*
Draws a cube
centered at (x,y,z)
with size (xLength,yLength,zLength)
rotated th degrees around (rx,ry,rz)
*/
void Scene::Cube (double x,double y,double z,
	   double xLength, double yLength, double zLength,
	   double th,double rx,double ry,double rz)
{
	glPushMatrix();
	// Transform 
	glTranslated(x,y,z);
	glRotated(th,rx,ry,rz);
   	glScaled(xLength/2,yLength/2,zLength/2);
	glBegin(GL_QUADS);
	// +x side
   	glNormal3f( 1, 0, 0);
	glVertex3d(1,1,1);
	glVertex3d(1,1,-1);
	glVertex3d(1,-1,-1);
	glVertex3d(1,-1,1);
	// +z side
  	glNormal3f( 0, 0, 1);
	glVertex3d(1,1,1);
	glVertex3d(-1,1,1);
	glVertex3d(-1,-1,1);
	glVertex3d(1,-1,1);
	// -x side
   	glNormal3f( -1, 0, 0);
	glVertex3d(-1,1,1);
	glVertex3d(-1,1,-1);
	glVertex3d(-1,-1,-1);
	glVertex3d(-1,-1,1);
	// -z side
  	glNormal3f( 0, 0, -1);
	glVertex3d(1,1,-1);
	glVertex3d(-1,1,-1);
	glVertex3d(-1,-1,-1);
	glVertex3d(1,-1,-1);
	// Top
   	glNormal3f( 0, 1, 0);
	glVertex3d(1,1,1);
	glVertex3d(-1,1,1);
	glVertex3d(-1,1,-1);
	glVertex3d(1,1,-1);
	// Bottom
	glNormal3f( 0, -1, 0);
	glVertex3d(1,-1,1);
	glVertex3d(-1,-1,1);
	glVertex3d(-1,-1,-1);
	glVertex3d(1,-1,-1);
	glEnd();
	glPopMatrix();
}

