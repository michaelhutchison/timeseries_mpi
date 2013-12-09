#ifndef CSCIx229
#define CSCIx229

//  Include GLEW if you need access to advanced features
#ifdef USEGLEW
#include <GL/glew.h>
#endif
//  Include GLUT and GL
#define GL_GLEXT_PROTOTYPES
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

//  Include standard headers
#include <stdio.h>
#include <stdarg.h>

//  Handy trigonometry macros 
#define PI 3.1415926
#define Cos(th) cos(PI/180*(th))
#define Sin(th) sin(PI/180*(th))

//  Functions provided by the library
void Print(const char* format , ...);
void Fatal(const char* format , ...);
void ErrCheck(const char* where);

#endif

