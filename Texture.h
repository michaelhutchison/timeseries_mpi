#ifndef TEXTURE_H
#define TEXTURE_H

#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include "graphicslib.h"

class Texture {
  public:
    Texture();
    ~Texture();
    unsigned int LoadTexBMP(const char* file);
    GLuint loadTexture(std::string source, std::string key);
    GLuint get(std::string key);
  private:
    void Reverse(void* x,const int n);
    std::map<std::string, GLuint>textureMap;
};

#endif
