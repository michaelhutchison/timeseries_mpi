#ifndef MOUSE_H
#define MOUSE_H

#include "graphicslib.h"
#include "View.h"

class Mouse {
  public:
    Mouse(View * v);
    ~Mouse();
    void motion(int x, int y);
    void button(int key,int status,int x,int y);
    int getx() {return x;}
    int gety() {return y;}
    int getMode() {return mode;}
  private:
    int mode;
    int x,y;
    int clickX,clickY;

    View * view;
};


#endif

