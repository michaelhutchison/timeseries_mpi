#include "vec3.h"



/* CONSTRUCTORS and DESTRUCTORS */
Vec3::Vec3() {
    x = v3_type();
    y = v3_type();
    z = v3_type();
}
Vec3::Vec3(v3_type x0, v3_type y0, v3_type z0) {
    x = x0;
    y = y0;
    z = z0;
}
Vec3::Vec3(v3_type xyz0) {
    x = xyz0;
    y = xyz0;
    z = xyz0;
} 

Vec3::Vec3( const Vec3& b){ 
    x = b.x;
    y = b.y;
    z = b.z;
}
Vec3::~Vec3() {

}

/* MEMBER FUNCTIONS */
Vec3 Vec3::cross(const Vec3& B) const{
    Vec3 v;
    v.x = y*B.z - z*B.y;
    v.y = z*B.x - x*B.z;
    v.z = x*B.y - y*B.x;
    return v;
}
Vec3::v3_type Vec3::dot(const Vec3& B) const{
    Vec3::v3_type m = x*B.x + y*B.y + z*B.z;
    return m;
}
Vec3::v3_type Vec3::length() const{
    return sqrt( x*x + y*y + z*z);
}
Vec3::v3_type Vec3::lengthSquared() const {
    return x*x + y*y + z*z;
}
Vec3::v3_type Vec3::angleDegrees(const Vec3& B) const{
    return (180.0 * angleRadians(B)) / PI;
}
Vec3::v3_type Vec3::angleRadians(const Vec3& B) const{
    return acos( dot(B) / (length() * B.length()));
}
Vec3 Vec3::normalize() const{
    Vec3::v3_type l = length();
    Vec3 v(x/l,y/l,z/l);
    return v;
}

/* MEMBER OPERATORS */
void Vec3::operator+=(const Vec3& b) {
    x+=b.x;
    y+=b.y;
    z+=b.z;    
}
void Vec3::operator-=(const Vec3& b){
    x-=b.x;
    y-=b.y;
    z-=b.z;
}
void Vec3::operator*=(const v3_type& c){
    x*=c;
    y*=c;
    z*=c;
}
void Vec3::operator/=(const v3_type& c){
    x/=c;
    y/=c;
    z/=c;
}

/* OPERATORS */
Vec3 operator+(Vec3 a, Vec3 b){
    Vec3 v;
    v.x = a.x + b.x;
    v.y = a.y + b.y;
    v.z = a.z + b.z;
    return v;
}
Vec3 operator-(Vec3 a, Vec3 b){
    Vec3 v;
    v.x = a.x - b.x;
    v.y = a.y - b.y;
    v.z = a.z - b.z;
    return v;
}
Vec3 operator*(Vec3 a, Vec3::v3_type c){
    Vec3 v;
    v.x = a.x * c;
    v.y = a.y * c;
    v.z = a.z * c;
    return v;
}
Vec3 operator*(Vec3::v3_type c, Vec3 a){
    return a * c;
}
Vec3 operator/(Vec3 a, Vec3::v3_type c){
    Vec3 v;
    v.x = a.x / c;
    v.y = a.y / c;
    v.z = a.z / c;
    return v;
}
std::ostream& operator<<(std::ostream& outs, const Vec3& a){
    outs << "(" << a.x << "," << a.y << "," << a.z << ")";
    return outs;
}

/* OPERATIONS */
Vec3 cross(const Vec3& a, const Vec3& b) {
    Vec3 v;
    
    
    v.x = (a.y*b.z) - (a.z*b.y);
    v.y = (a.z*b.x) - (a.x*b.z);
    v.z = (a.x*b.y) - (a.y*b.x);
    return v;
}
Vec3::v3_type dot(const Vec3& a,const Vec3& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
Vec3::v3_type angleDegrees(const Vec3& a, const Vec3& b){
    return a.angleDegrees(b);
}
Vec3::v3_type angleRadians(const Vec3& a, const Vec3& b){
    return a.angleRadians(b);
}
Vec3 normalize(const Vec3& a) {
    return a/a.length();
}

