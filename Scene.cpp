#include "Scene.h"


Scene::Scene() : model(1000) {
    axesOff();
    W = 600;
    H = 600;
    backgroundColor[0] = 1.0;
    backgroundColor[1] = 1.0;
    backgroundColor[2] = 1.0;
    backgroundColor[3] = 1.0;

    // initialize colors
    colors[ 0].set(1.0, 0.0, 0.0);
    colors[ 1].set(0.0, 1.0, 0.0);
    colors[ 2].set(0.0, 0.0, 1.0);
    colors[ 3].set(1.0, 1.0, 0.0);
    colors[ 4].set(1.0, 0.0, 1.0);
    colors[ 5].set(0.0, 1.0, 1.0);
    colors[ 6].set(1.0, 0.5, 0.0);
    colors[ 7].set(1.0, 0.0, 0.5);
    colors[ 8].set(0.5, 1.0, 0.0);
    colors[ 9].set(0.0, 1.0, 0.5);
    colors[10].set(0.5, 0.0, 1.0);
    colors[11].set(0.0, 0.5, 1.0);
    colors[12].set(1.0, 1.0, 0.5);
    colors[13].set(1.0, 0.5, 1.0);
    colors[14].set(0.5, 1.0, 1.0);
    colors[15].set(1.0, 0.1, 0.3);

    // Initialize the light
    light.setPosition(60,60,60);
}
Scene::~Scene() {

}
void Scene::initialize(char * file) {
    strncpy(filename, file, 90);
    loadTextures();
    reader.read(filename);
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
    if (drawAxesOn) drawAxes(-55,-55,-55);

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
    // Draw box bounding the world
    WireframeCube(0,0,0,
                  100,100,100,
                  0,0,0,0);

    // Draw recorded scene
    double radius = 1.0;
    FrameState * frame = reader.get_next_frame();
    if (frame != 0) {
        unsigned long nObjects = frame->get_num_objects();
        for (unsigned long i=0; i<nObjects; i++) {
            Object * o = frame->get_object(i);
            if (o != 0) {
                // Set color based on id
                int origin = o->getID() >> 28;
 	            glColor3f(colors[origin].x, colors[origin].y, colors[origin].z);


                // Draw an object from the animation record file
                //Cube(o->x(), o->y(), o->z(),
                //           0.4,0.4,0.4,
                //           o->a_r(),o->x_r(),o->y_r(),o->z_r());
                Sphere(o->x(), o->y(), o->z(), radius);
            }
        }
    }
}
/*
Draws X Y Z axes
*/
void Scene::drawAxes(double x, double y, double z)
{
    glPushMatrix();
	glTranslated(x,y,z);
	double len=5;
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
	glPopMatrix();
}

/*
Draws a wirefram cube
centered at (x,y,z)
with size (xLength,yLength,zLength)
rotated th degrees around (rx,ry,rz)
*/
void Scene::WireframeCube (double x,double y,double z,
	   double xLength, double yLength, double zLength,
	   double th,double rx,double ry,double rz)
{
	glPushMatrix();
	// Transform 
	glTranslated(x,y,z);
	glRotated(th,rx,ry,rz);
   	glScaled(xLength/2,yLength/2,zLength/2);
	// +x side and -x side and one edge
	glBegin(GL_LINE_STRIP);
	// +x side and one edge
    glVertex3d(1,1,1);
	glVertex3d(1,1,-1);
	glVertex3d(1,-1,-1);
	glVertex3d(1,-1,1);
	glVertex3d(1,1,1);
    // -x side
	glVertex3d(-1,1,1);
	glVertex3d(-1,1,-1);
	glVertex3d(-1,-1,-1);
	glVertex3d(-1,-1,1);
	glVertex3d(-1,1,1);
    glEnd();
    // remaining edges parallel to x-axis
    glBegin(GL_LINES);
	glVertex3d( 1, 1,-1);
	glVertex3d(-1, 1,-1);
	glVertex3d( 1,-1,-1);
	glVertex3d(-1,-1,-1);
	glVertex3d( 1,-1, 1);
	glVertex3d(-1,-1, 1);
    glEnd();
	glPopMatrix();
}
void Scene::Sphere(double x, double y, double z, double radius) {
	double inc=45;
	glPushMatrix();
	glTranslated(x,y,z);
	glScaled(radius,radius,radius);
   	// Begin drawing
	glBegin(GL_QUAD_STRIP);	
	for (double i=0;i<=180;i+=inc)
	{
		for (double j=0;j<=360;j+=inc)
		{
   			glNormal3f( Cos(j+0  )*Sin(i+0  ), Cos(i+0  ), Sin(j+0  )*Sin(i+0  ));
			glTexCoord2f( (j+0  )/360,(i+0  )/180);
			glVertex3f( Cos(j+0  )*Sin(i+0  ), Cos(i+0  ), Sin(j+0  )*Sin(i+0  ));
			glNormal3f( Cos(j+0  )*Sin(i+inc), Cos(i+inc), Sin(j+0  )*Sin(i+inc));
			glTexCoord2f( (j+0  )/360,(i+inc)/180);
			glVertex3f( Cos(j+0  )*Sin(i+inc), Cos(i+inc), Sin(j+0  )*Sin(i+inc));
		}
	}
	glEnd();
	// Finished drawing	
    glUseProgram(0);
	glPopMatrix();
   	glDisable(GL_TEXTURE_2D);
}


