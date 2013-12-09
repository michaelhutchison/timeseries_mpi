#ifndef VEC3_H
#define VEC3_H

#ifndef PI
#define PI 3.14159265
#endif

#include <math.h>
#include <iostream>

class Vec3 {
  public:
    typedef double v3_type;
    /* PUBLIC DATA */
    v3_type x;
    v3_type y;
    v3_type z;
    
    /* CONSTRUCTORS AND DESTRUCTOR */
    Vec3();
    Vec3(v3_type x0, v3_type y0, v3_type z0);
    Vec3(v3_type xyz0);
    Vec3(const Vec3& b); // copy constructor
    ~Vec3();

    /* PUBLIC METHODS */
    void set(v3_type a, v3_type b, v3_type c) {x=a; y=b; z=c;}
    void set(v3_type a, v3_type b) {x=a; y=b;}
    Vec3 cross(const Vec3& B) const;
    v3_type dot(const Vec3& B) const;
    v3_type length() const;
    v3_type lengthSquared() const;
    v3_type angleDegrees(const Vec3& B) const;
    v3_type angleRadians(const Vec3& B) const;
    Vec3 normalize() const;
  
    void operator+=(const Vec3& b);
    void operator-=(const Vec3& b);
    void operator*=(const v3_type& b);
    void operator/=(const v3_type& b);

    /* ALIASES */
    v3_type r()     {return x;}
    v3_type g()     {return y;}
    v3_type b()     {return z;}
    v3_type roll()  {return x;}
    v3_type pitch() {return y;}
    v3_type yaw()   {return z;}
  private:

};

/* OPERATORS */
Vec3 operator+(Vec3 a, Vec3 b);
Vec3 operator-(Vec3 a, Vec3 b);
Vec3 operator*(Vec3 a, Vec3::v3_type c);
Vec3 operator*(Vec3::v3_type c, Vec3 a);
Vec3 operator/(Vec3 a, Vec3::v3_type c);
std::ostream& operator<<(std::ostream& outs, const Vec3& a);

/* OPERATIONS */
Vec3 cross(const Vec3& a, const Vec3& b);
Vec3::v3_type dot(const Vec3& a,const Vec3& b);
Vec3::v3_type angleDegrees(const Vec3& a,const Vec3& b);
Vec3::v3_type angleRadians(const Vec3& a,const Vec3& b);
Vec3 normalize(const Vec3& a);

#endif
