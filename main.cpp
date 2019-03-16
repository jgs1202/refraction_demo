// [header]
// A very basic raytracer example.
// [/header]
// [compile]
// c++ -o raytracer -O3 -Wall raytracer.cpp
// [/compile]
// [ignore]
// Copyright (C) 2012  www.scratchapixel.com
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// [/ignore]
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <cfenv>
#include <stdlib.h>
#include <math.h>
#include <fenv.h>

#include <GLUT/glut.h>
#include <openGL/gl.h>
#include <openGL/glu.h>


#if defined __linux__ || defined __APPLE__
// "Compiled for Linux
#else
// Windows doesn't define these values by default, Linux does
#define M_PI 3.141592653589793
#define INFINITY 1e8
#endif

// each length must be the same because of the program
#define WIDTH 100
#define HEIGHT 100
#define DEPTH 100
//#define WIDTH 5
//#define HEIGHT 5
//#define DEPTH 5
#define STEP 1
#define MAX_RAY_DEPTH 5

template<typename T>
class Vec3
{
public:
    T x, y, z;
    Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
    Vec3(T xx) : x(xx), y(xx), z(xx) {}
    Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
    Vec3& normalize()
    {
        T nor2 = length2();
        if (nor2 > 0) {
            T invNor = 1 / sqrt(nor2);
            x *= invNor, y *= invNor, z *= invNor;
        }
        return *this;
    }
    Vec3<T> operator * (const T &f) const { return Vec3<T>(x * f, y * f, z * f); }
    Vec3<T> operator / (const T &f) const { return Vec3<T>(x / f, y / f, z / f); }
    Vec3<T> operator * (const Vec3<T> &v) const { return Vec3<T>(x * v.x, y * v.y, z * v.z); }
    T dot(const Vec3<T> &v) const { return x * v.x + y * v.y + z * v.z; }
    Vec3<T> operator - (const Vec3<T> &v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); }
    Vec3<T> operator + (const Vec3<T> &v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); }
    Vec3<T>& operator += (const Vec3<T> &v) { x += v.x, y += v.y, z += v.z; return *this; }
    Vec3<T>& operator *= (const Vec3<T> &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
    Vec3<T> operator - () const { return Vec3<T>(-x, -y, -z); }
    T length2() const { return x * x + y * y + z * z; }
    T length() const { return sqrt(length2()); }
    friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v)
    {
        os << "[" << v.x << " " << v.y << " " << v.z << "]";
        return os;
    }
};

typedef Vec3<float> Vec3f;

typedef Vec3<int> Vec3i;

class Sphere
{
public:
    Vec3f center;                           /// position of the sphere
    float radius, radius2;                  /// sphere radius and radius^2
    Vec3f surfaceColor, emissionColor;      /// surface color and emission (light)
    float transparency, reflection;         /// surface transparency and reflectivity
    Sphere(
            const Vec3f &c,
            const float &r,
            const Vec3f &sc,
            const float &refl = 0,
            const float &transp = 0,
            const Vec3f &ec = 0) :
            center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec),
            transparency(transp), reflection(refl)
    { /* empty */ }
    //[comment]
    // Compute a ray-sphere intersection using the geometric solution
    //[/comment]
    bool intersect(const Vec3f &rayorig, const Vec3f &raydir, float &t0, float &t1) const
    {
        Vec3f l = center - rayorig;
        float tca = l.dot(raydir);
        if (tca < 0) return false;
        float d2 = l.dot(l) - tca * tca;
        if (d2 > radius2) return false;
        float thc = sqrt(radius2 - d2);
        t0 = tca - thc;
        t1 = tca + thc;

        return true;
    }
};

struct object
{
public:
    Vec3f coordinate;
    Vec3f color;
    float reflectivity;
    float opacity;
    float refractivity;
    Vec3f emissionColor;
};

struct light
{
    Vec3f direction;
    Vec3f color;
    float intensity;
};


// declare variables here as global
object field[WIDTH * HEIGHT * DEPTH];
light fromLightSource[WIDTH * HEIGHT * DEPTH];
Vec3f grad[WIDTH * HEIGHT * DEPTH];
Vec3f background = (1);

float roundInt (float x){
    float decimal = x - std::floor(x);
    if (decimal < 0.5) {
        return std::floor(x);
    } else {
        return std::ceil(x);
    }
}

float vectorDistance (Vec3f vec1, Vec3f vec2){
    return sqrt(pow((vec1.x - vec2.x), 2) + pow((vec1.y - vec2.y), 2) + pow((vec1.z - vec2.z), 2));
}

float slope2 (float p1, float p2){
    return (p2 - p1) / STEP;
}
float slope3 (float p1, float p2, float p3){
    return (p2 - p1) / STEP / 2 + (p3 - p2) / STEP / 2;
}

void gradRefraction(object *field, Vec3f grad[]){
//    Vec3f grad[WIDTH * HEIGHT * DEPTH];
    for (unsigned int i=0; i<WIDTH; ++i){
        for (unsigned int j=0; j<HEIGHT; ++j){
            for (unsigned int k=0; k<DEPTH; ++k){
                float x1 = 0., x2 = 0., y1 = 0., y2 = 0., z1 = 0., z2 = 0.;
                bool xx1 = false, xx2 = false, yy1 = false, yy2 = false, zz1 = false, zz2 = false;
                float origin = field[i + j * WIDTH + k * WIDTH * HEIGHT].refractivity;
                if (i != 0) {
                    xx1 = true;
                    x1 = field[i - 1 + j * WIDTH + k * WIDTH * HEIGHT].refractivity;
                }
                if (i != WIDTH - 1) {
                    xx2 = true;
                    x2 = field[i + 1 + j * WIDTH + k * WIDTH * HEIGHT].refractivity;
                } if (j != 0){
                    yy1 = true;
                    y1 = field[i + (j - 1) * WIDTH + k * WIDTH * HEIGHT].refractivity;
                } if (j != HEIGHT - 1) {
                    yy2 = true;
                    y2 = field[i + (j + 1) * WIDTH + k * WIDTH * HEIGHT].refractivity;
                } if (k != 0){
                    zz1 = true;
                    z1 = field[i + j * WIDTH + (k - 1) * WIDTH * HEIGHT].refractivity;
                } if (k != DEPTH - 1) {
                    zz2 = true;
                    z2 = field[i + j * WIDTH + (k + 1) * WIDTH * HEIGHT].refractivity;
                }
                grad[i + j * WIDTH + k * WIDTH * HEIGHT].x = (xx1 * xx2)? slope3(x1, origin, x2) : ((xx1)? slope2(x1, origin) : slope2(origin, x2));
                grad[i + j * WIDTH + k * WIDTH * HEIGHT].y = (yy1 * yy2)? slope3(y1, origin, y2) : ((yy1)? slope2(y1, origin) : slope2(origin, y2));
                grad[i + j * WIDTH + k * WIDTH * HEIGHT].z = (zz1 * zz2)? slope3(z1, origin, z2) : ((zz1)? slope2(z1, origin) : slope2(origin, z2));

//                std::cout << xx1 << ", " << xx2 << ",  " << ((xx1 * xx2)? (2) : ((xx1)? (1) : (0))) << std::endl;
            }
        }
    }
}

// position, radius, surface color, reflectivity, transparency, refractivity, emission color
void inSphere(object field[],  Vec3f center, float radius, Vec3f surfaceColor, float reflectivity, float opacity, float refractivity, Vec3f emissionColor){
    for (unsigned int i=0; i < WIDTH; ++i){
        for (unsigned int j=0; j<HEIGHT; ++j){
            for (unsigned int k=0; k<DEPTH; ++k){
                if (std::abs(i - center.x) > radius){
                    continue;
                } else {
                } if (std::abs(j - center.y) > radius){
                    continue;
                } if (std::abs(k - center.z) > radius){
                    continue;
                } if (((i - center.x) * (i - center.x) + (j - center.y) * (j - center.y) + (k - center.z) * (k - center.z)) > radius * radius){
                    continue;
                }
                field[i + j * WIDTH + k * WIDTH * HEIGHT].coordinate = center;
                field[i + j * WIDTH + k * WIDTH * HEIGHT].color = surfaceColor;
                field[i + j * WIDTH + k * WIDTH * HEIGHT].reflectivity = reflectivity;
                field[i + j * WIDTH + k * WIDTH * HEIGHT].opacity = opacity;
                field[i + j * WIDTH + k * WIDTH * HEIGHT].emissionColor = emissionColor;
                field[i + j * WIDTH + k * WIDTH * HEIGHT].refractivity = refractivity;
            }
        }
    }
}

unsigned int getCoordinate(Vec3f coordinate){
    return coordinate.x + coordinate.y * WIDTH + coordinate.z * WIDTH * WIDTH;
}

Vec3f vecNormalize(Vec3f vec){
    Vec3f sign = (vec.x / std::abs(vec.x), vec.y / std::abs(vec.y), vec.z / std::abs(vec.z));
    float square = vec.x * vec.x + vec.y * vec.y + vec.z + vec.z;
    return (vec / std::sqrt(square)) * sign;
}

void inPlain (Vec3f p1, Vec3f p2, Vec3f p3, Vec3f p4, Vec3f color, float reflectivity, float opacity, float refractivity){
    unsigned int dis1 = std::floor(vectorDistance(p1, p2)) + 1;
    unsigned int dis2 = std::floor(vectorDistance(p1, p3)) + 1;
    Vec3f vec1 = vecNormalize(p2 - p1);
    Vec3f vec2 = vecNormalize(p3 - p1);
    Vec3f current;
    for (unsigned int x=0; x < dis1; ++x){
        for (unsigned y=0; y < dis2; ++y){
            current = p1 + vec1 * x + vec2 * y;
            current = Vec3f(roundInt(current.x), roundInt(current.y), roundInt(current.z));
            field[getCoordinate(current)].color = color;
            field[getCoordinate(current)].reflectivity = reflectivity;
            field[getCoordinate(current)].opacity = opacity;
            field[getCoordinate(current)].refractivity = refractivity;

        }
    }
}

Vec3f getVector(unsigned int coordinate){
    unsigned int x, y, z;
    Vec3f vector;
    z = (int) coordinate / (WIDTH * WIDTH);
    y = (int) (coordinate - z * WIDTH * WIDTH) / WIDTH;
    x = (int) (coordinate - y * WIDTH  - z * WIDTH * WIDTH);
    vector = Vec3f (x, y, z);
    return vector;
}

Vec3f identityVector(unsigned int direction){
    Vec3f initial;
    switch(direction){
        case 0:
            initial = Vec3f(1., 0., 0.);
            break;
        case 1:
            initial = Vec3f(0., 1., 0.);
            break;
        case 2:
            initial = Vec3f(0., 0., 1.);
            break;
    }
    return initial;
}

Vec3f absVec3(Vec3f vec){
    return Vec3f(std::abs(vec.x), std::abs(vec.y), std::abs(vec.z));
}

float inpro(Vec3f vec1, Vec3f vec2){
    return vec1.x * vec2.x + vec1.y + vec2.y + vec1.z + vec2.z;
}

Vec3f checkMinus(Vec3f vec){
    Vec3f vec2 = vec;
    if (vec.x < 0){
        vec2.x = 0;
    } if (vec.y < 0) {
        vec2.y = 0;
    } if (vec.z < 0) {
        vec2.z = 0;
    }
    return vec2;
}

Vec3f bilinearInterpolation(Vec3f origin, Vec3f p1, Vec3f p2, Vec3f p3, Vec3f p4, light light[]){
//    std::cout << p1 << p2 << p3 << p4 << std::endl;
    Vec3f x1, x2, y1, y2;
    Vec3f xx1, xx2;
    x1 = p1;
    float dis1 = vectorDistance(p1, p2);
    float dis2 = vectorDistance(p1, p3);
    float dis3 = vectorDistance(p1, p4);
    if ((dis1 > dis2) || (dis1 > dis3) ) {
        x2 = p3;
        y1 = p4;
        y2 = p2;
    } else if ((dis2 > dis1) || (dis2 > dis3)){
        y2 = p3;
        x2 = p2;
        y1 = p4;
    } else {
        y2 = p4;
        x2 = p2;
        y1 = p3;
    }
    // get direction vector like (0, 0, 1)
    Vec3f dic1 = x1 - x2;
    Vec3f dic2 = x1 - y1;
    Vec3f ratio1 = absVec3(origin * dic1 - x1 * dic1) / STEP;
    Vec3f ratio2 = absVec3(origin * dic2 - x1 * dic2 ) / STEP;

    // 座標がマイナスなら保管
    x1 = checkMinus(x1);
    x2 = checkMinus(x2);
    y1 = checkMinus(y1);
    y2 = checkMinus(y2);
    xx1 = light[getCoordinate(x1)].direction * (1 - inpro(Vec3f(1, 1, 1), ratio1)) + light[getCoordinate(x2)].direction + inpro(Vec3f(1, 1, 1), ratio1);
    xx2 = light[getCoordinate(y1)].direction * (1 - inpro(Vec3f(1, 1, 1), ratio1)) + light[getCoordinate(y2)].direction + inpro(Vec3f(1, 1, 1), ratio1);
    xx1 = xx1.normalize();
    xx2 = xx2.normalize();
    return xx1 * (1 - inpro(Vec3f(1, 1, 1), ratio2)) + xx2 * inpro(Vec3f(1, 1, 1), ratio2);
}

bool checkInteger(Vec3f vec){
    if ((std::ceil(vec.x) == std::floor(vec.x)) && (std::ceil(vec.y) == std::floor(vec.y)) && (std::ceil(vec.z) == std::floor(vec.z)) ){
        return true;
    } else {
        return false;
    }
}


//
//Vec3f nextDirection(light light[], unsigned int coordinate, unsigned int previousCoo, Vec3f grad[], unsigned int direction) {
//    Vec3f previousDirection = light[previousCoo].direction + grad[previousCoo] * STEP * (-1);
//    Vec3f backPoint, p1, p2, p3, p4;
//    float coeff;
//
//    if (direction == 0) {
//        coeff = 1 / previousDirection.x;
//        backPoint = getVector(coordinate) - previousDirection * coeff;
//        backPoint.y = (unsigned int) backPoint.y;
//        backPoint.z = (unsigned int) backPoint.z;
//        p1 = getVector(coordinate) - backPoint;
//        p2 = p1 - identityVector(1);
//        p3 = p2 - identityVector(2);
//        p4 = p1 - identityVector(2);
//    } else if (direction == 1) {
//        coeff = 1 / previousDirection.y;
//        backPoint = getVector(coordinate) - previousDirection * coeff;
//        backPoint.z = (unsigned int) backPoint.z;
//        backPoint.x = (unsigned int) backPoint.x;
//        p1 = getVector(coordinate) - backPoint;
//        p2 = p1 - identityVector(0);
//        p3 = p2 - identityVector(2);
//        p4 = p1 - identityVector(2);
//    } else if (direction == 2) {
//        coeff = 1 / previousDirection.z;
//        backPoint = getVector(coordinate) - previousDirection * coeff;
//        backPoint.x = (unsigned int) backPoint.x;
//        backPoint.y = (unsigned int) backPoint.y;
//        p1 = getVector(coordinate) - backPoint;
//        p2 = p1 - identityVector(0);
//        p3 = p2 - identityVector(1);
//        p4 = p1 - identityVector(1);
//    }
//    if (checkInteger) {
//        return previousDirection;
//    } else {
//        return bilinearInterpolation(getVector(coordinate) - previousDirection * coeff, p1, p2, p3, p4, light);
//    }
//}

//void calcLightSource(object *field, light light[], Vec3f grad[], unsigned int direction, Vec3f color, float intensity){
//    unsigned int coordinate;
//    unsigned int previousCoo;
//    for (unsigned int i=0; i<WIDTH; ++i){
//       for (unsigned int j=0; j<HEIGHT; ++j){
//           for (unsigned int k=0; k<DEPTH; ++k){
//               coordinate = i * (int)pow(WIDTH, direction) + j * (int)pow(WIDTH, (direction==0)?2:(direction - 1)) + k * (int)pow(WIDTH, (direction<2)?(direction + 1):(0));
//               if (i == 0){
//                   Vec3f initial = identityVector(direction);
//                   light[coordinate].color = color;
//                   light[coordinate].intensity = intensity;
//                   light[coordinate].direction = initial.normalize();
//               } else {
//                   previousCoo = coordinate - (int)pow(WIDTH, direction);
//                   light[coordinate].color = light[previousCoo].color * field[previousCoo].color * (1 - field[previousCoo].opacity) * light[previousCoo].intensity;
//                   light[coordinate].direction = nextDirection(light, coordinate, previousCoo, grad, direction).normalize();
//                   light[coordinate].intensity = light[previousCoo].intensity * (1 - field[previousCoo].opacity);
//               }
//           }
//       }
//    }
//}

Vec3f nextDirection(light light[], unsigned int coordinate, unsigned int previousCoo, Vec3f grad[], Vec3f direction,
                     unsigned int maxDir) {
    Vec3f previousDirection = light[previousCoo].direction + grad[previousCoo] * STEP * (-1);
    Vec3f backPoint, p1, p2, p3, p4;
    float coeff;

    if (maxDir == 0) {
        coeff = 1 / previousDirection.x;
        backPoint = getVector(coordinate) - previousDirection * coeff;
        backPoint.y = (unsigned int) backPoint.y;
        backPoint.z = (unsigned int) backPoint.z;
        p1 = getVector(coordinate) - backPoint;
        p2 = p1 - identityVector(1);
        p3 = p2 - identityVector(2);
        p4 = p1 - identityVector(2);
    } else if (maxDir == 1) {
        coeff = 1 / previousDirection.y;
        backPoint = getVector(coordinate) - previousDirection * coeff;
        backPoint.z = (unsigned int) backPoint.z;
        backPoint.x = (unsigned int) backPoint.x;
        p1 = getVector(coordinate) - backPoint;
        p2 = p1 - identityVector(0);
        p3 = p2 - identityVector(2);
        p4 = p1 - identityVector(2);
    } else if (maxDir == 2) {
        coeff = 1 / previousDirection.z;
        backPoint = getVector(coordinate) - previousDirection * coeff;
        backPoint.x = (unsigned int) backPoint.x;
        backPoint.y = (unsigned int) backPoint.y;
        p1 = getVector(coordinate) - backPoint;
        p2 = p1 - identityVector(0);
        p3 = p2 - identityVector(1);
        p4 = p1 - identityVector(1);
    }
    if (checkInteger) {
        return previousDirection;
    } else {
        return bilinearInterpolation(getVector(coordinate) - previousDirection * coeff, p1, p2, p3, p4, light);
    }
}

void calcLightSource(object *field, light light[], Vec3f grad[], Vec3f direction, Vec3f color, float intensity){
    unsigned int coordinate;
    unsigned int previousCoo;
    unsigned int maxDir;
    if ((direction.x >= direction.y) && (direction.x >= direction.z)){
        maxDir = 0;
    } else if ((direction.y >= direction.x) && (direction.y >= direction.z)){
        maxDir = 1;
    } else {
        maxDir = 2;
    }
    for (unsigned int i=0; i<WIDTH; ++i){
        for (unsigned int j=0; j<HEIGHT; ++j){
            for (unsigned int k=0; k<DEPTH; ++k){
                coordinate = i * (int)pow(WIDTH, maxDir) + j * (int)pow(WIDTH, (maxDir==0)?2:(maxDir - 1)) + k * (int)pow(WIDTH, (maxDir<2)?(maxDir + 1):(0));
                if (i == 0){
                    light[coordinate].color = color;
                    light[coordinate].intensity = intensity;
                    light[coordinate].direction = vecNormalize(vecNormalize(direction));
                } else {
                    previousCoo = coordinate - (int)pow(WIDTH, maxDir);
                    light[coordinate].color = light[previousCoo].color * field[previousCoo].color * (1 - field[previousCoo].opacity) * light[previousCoo].intensity;
                    light[coordinate].direction = vecNormalize(nextDirection(light, coordinate, previousCoo, grad, direction, maxDir));
                    light[coordinate].intensity = light[previousCoo].intensity * (1 - field[previousCoo].opacity);
                }
            }
        }
    }
}

float mix(const float &a, const float &b, const float &mix)
{
    return b * mix + a * (1 - mix);
}

//[comment]
// This is the main trace function. It takes a ray as argument (defined by its origin
// and direction). We test if this ray intersects any of the geometry in the scene.
// If the ray intersects an object, we compute the intersection point, the normal
// at the intersection point, and shade this point using this information.
// Shading depends on the surface property (is it transparent, reflective, diffuse).
// The function returns a color for the ray. If the ray intersects an object that
// is the color of the object at the intersection point, otherwise it returns
// the background color.
//[/comment]

void bilinearTrace(float xx, float yy, unsigned int z, Vec3f& color, Vec3f& medium, Vec3f& reflectivity, float& opacity, Vec3f& gradient, float& intensity, bool& flagEnd){
    bool flagInt = false;
    unsigned int p;
    if ((std::ceil(xx) == std::floor(xx))){
        if ((std::ceil(yy) == std::floor(yy))) {
            flagInt = true;
        }
    }
    if ((xx < 0) || (yy < 0) || (xx >= WIDTH) || ((yy >= HEIGHT))){
//        std::cout << "case 1" << std::endl;
        xx = 0, yy = 0;

        flagEnd = true;
        p = (int)std::floor(xx) + (int)std::floor(yy) * WIDTH + z * WIDTH * HEIGHT;
        color = fromLightSource[p].color;
        medium = field[p].color;
        reflectivity = field[p].reflectivity;
        opacity = field[p].opacity;
        gradient = grad[p];
        intensity = fromLightSource[p].intensity;
    } else if (flagInt) {
//        std::cout << "case 2 " << xx << " " << yy << " " << z << std::endl;
        p = (int)xx + (int)yy * WIDTH + z * WIDTH * HEIGHT;
        color = fromLightSource[p].color;
        medium = field[p].color;
        reflectivity = field[p].reflectivity;
        opacity = field[p].opacity;
        gradient = grad[p];
        intensity = fromLightSource[p].intensity;
    } else {
//        std::cout << "case 3" << z << std::endl;
        Vec3f x1, x2, y1, y2;
        x1 = Vec3f((int)std::floor(xx), (int)std::floor(yy), z);
        x2 = Vec3f((int)std::floor(xx), (int)std::floor(yy) + 1, z);
        y1 = Vec3f((int)std::floor(xx) + 1, (int)std::floor(yy), z);
        y2 = Vec3f((int)std::floor(xx) + 1, (int)std::floor(yy) + 1, z);

        float ratiox = std::abs(xx - std::floor(xx)) / STEP;
        float ratioy = std::abs(yy - std::floor(yy)) / STEP;

        color = (fromLightSource[getCoordinate(x1)].color * (1 - ratiox) + fromLightSource[getCoordinate(y1)].color * ratiox) * (1 - ratioy) + (fromLightSource[getCoordinate(x2)].color * (1 - ratiox) + fromLightSource[getCoordinate(y2)].color * ratiox) * ratioy;
        medium = (field[getCoordinate(x1)].color * (1 - ratiox) + field[getCoordinate(y1)].color * ratiox) * (1 - ratioy) + (field[getCoordinate(x2)].color * (1 - ratiox) + field[getCoordinate(y2)].color * ratiox) * ratioy;
        reflectivity = (field[getCoordinate(x1)].reflectivity * (1 - ratiox) + field[getCoordinate(y1)].reflectivity * ratiox) * (1 - ratioy) + (field[getCoordinate(x2)].reflectivity * (1 - ratiox) + field[getCoordinate(y2)].reflectivity * ratiox) * ratioy;
        opacity = (field[getCoordinate(x1)].opacity * (1 - ratiox) + field[getCoordinate(y1)].opacity * ratiox) * (1 - ratioy) + (field[getCoordinate(x2)].opacity * (1 - ratiox) + field[getCoordinate(y2)].opacity * ratiox) * ratioy;
        gradient = (grad[getCoordinate(x1)] * (1 - ratiox) + grad[getCoordinate(y1)] * ratiox) * (1 - ratioy) + (grad[getCoordinate(x2)] * (1 - ratiox) + grad[getCoordinate(y2)] * ratiox) * ratioy;
        intensity = (fromLightSource[getCoordinate(x1)].intensity * (1 - ratiox) + fromLightSource[getCoordinate(y1)].intensity * ratiox) * (1 - ratioy) + (fromLightSource[getCoordinate(x2)].intensity * (1 - ratiox) + fromLightSource[getCoordinate(y2)].intensity * ratiox) * ratioy;
    }
}

Vec3f trace(Vec3f angle, unsigned int x, unsigned int y) {
    float xx = x, yy = y;
    unsigned int z = 0;
    Vec3f color = Vec3f(0., 0., 0.), medium = Vec3f(1.0, 1.0, 1.0);
    float opacity = 0;
    Vec3f zero = 0.;

    while (z < DEPTH){
//        update color, ld,
        Vec3f newColor, newMedium, newReflectivity, newGradient;

        float newIntensity, newOpacity;
        bool flagEnd = false;
        bilinearTrace(xx, yy, z, newColor, newMedium, newReflectivity, newOpacity, newGradient, newIntensity, flagEnd);
        if (flagEnd){
            break;
        }
        medium *= newMedium;
        opacity += (1 - opacity) * newOpacity;
        color += medium * (newColor * opacity) * (1 - opacity);
        if (std::abs(opacity - 1.0) < 0.1){
            break;
        }
        angle += newGradient * STEP;
        angle = angle * (1 / angle.z);
        xx += angle.x;
        yy += angle.y;
        z++;
    }
    return color;
}

// view point is always (0, 0, -distance)
Vec3f *render(object field[], light light[], float distance) {
    Vec3f *image = new Vec3f[WIDTH * HEIGHT], *pixel = image;
//    float invWidth = 1 / float(WIDTH), invHeight = 1 / float(HEIGHT);
//    float fov = 30, aspectratio = WIDTH / float(HEIGHT);
//    float angle = tan(M_PI * 0.5 * fov / 180.);
    Vec3f angle;
    // Trace rays
    for (int y = 0; y < HEIGHT; ++y) {
        for (int x = 0; x < WIDTH; ++x) {
            angle = Vec3f(x - WIDTH/2, y - HEIGHT/2, distance);
            angle *= 1 / angle.z;
            pixel[x + y * HEIGHT] = trace(angle, x, y);
        }
    }
//    for (int i=0; i<WIDTH * HEIGHT ; ++i){
//        std::cout << image[i] << std::endl;
//    }
    return pixel;
}

void myInit (char *progname)
{
    float aspect = (float) WIDTH / (float) HEIGHT;
    glutInitWindowSize( WIDTH, HEIGHT );
    glutInitDisplayMode( GLUT_SINGLE | GLUT_RGBA | GLUT_DEPTH);
//    glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutCreateWindow(progname);
    glClearColor (0.0, 0.0, 0.0, 1.0);
//    座標反転
    glPixelZoom(1.0f, -1.0f);
    glLoadIdentity();
}

void display( void )
{
//  座標調整
    glRasterPos2i(-1 , 1);
    glClear(GL_COLOR_BUFFER_BIT);

    std::cout << "display" << std::endl;

//    Vec3f *image = new Vec3f[WIDTH * HEIGHT], *pixel = image;
//    for (unsigned int i=0; i< 100 * 100; ++i){
//        pixel[i] = Vec3f(1, 1, 1);
//    }

//    glDrawPixels(WIDTH, HEIGHT, GL_RGB, GL_FLOAT, image);
    std::cout << "draw" << std::endl;
    glDrawPixels(WIDTH, HEIGHT, GL_RGB, GL_FLOAT, render(field, fromLightSource, 30));
    glFlush();
    std::cout << "display end" << std::endl;
}

void myReshape(int w, int h)
{
//    float aspect = (float) WIDTH / (float) HEIGHT;
    glViewport(0, 0, w, h);
//    glMatrixMode(GL_PROJECTION);
//    glLoadIdentity();
//    gluPerspective(30.0, aspect, 1.0, 50.0);
//    glMatrixMode(GL_MODELVIEW);
}

//[comment]
// In the main function, we will create the scene which is composed of 5 spheres
// and 1 light (which is also a sphere). Then, once the scene description is complete
// we render that scene, by calling the render() function.
//[/comment]

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    myInit(argv[0]);

    for (unsigned int i=0; i < WIDTH * HEIGHT * DEPTH; ++i){
        field[i] = {0};
        field[i].color = Vec3f(1.0, 1.0, 1.0);
    }


//    field, center, radius, color, reflection, opacity, refraction, emissionColor
    inSphere(field,  Vec3f (roundInt(WIDTH / 2), roundInt(WIDTH / 2), roundInt(WIDTH / 2)), roundInt(WIDTH / 3), Vec3f (1., 0., 0.), 0., 0.5, 0.5, Vec3f(0., 0., 0.));
    inPlain(Vec3f(0, HEIGHT, 0), Vec3f(WIDTH, HEIGHT, 0), Vec3f(0, HEIGHT, DEPTH), Vec3f(WIDTH, HEIGHT, DEPTH), Vec3f(1., 1., 1.), 0.3, 1., 0.);
    std::cout << field[getCoordinate(Vec3f(30, HEIGHT, 30))].color;
    gradRefraction(field, grad);

//    for (int i=0; i<WIDTH * HEIGHT; ++i){
//        std::cout << field[i + 2 * WIDTH * HEIGHT].color << std::endl;
//    }
    std::cout << std::endl << std::endl;

//  field, light, direction, color, intensity
    calcLightSource(field, fromLightSource, grad, Vec3f (1, 1, 0), Vec3f (1,1,1), 1.0);

//    unsigned int position = 14021;
//    for (int position=0; position<WIDTH*HEIGHT*DEPTH; ++position){
//        std::cout << getVector(position) << fromLightSource[position].color << ", " << fromLightSource[position].direction << ", " << fromLightSource[position].intensity << std::endl;
//    }
//
    glutDisplayFunc(display);
    glutReshapeFunc (myReshape);
    glutMainLoop();

    return 0;
}
