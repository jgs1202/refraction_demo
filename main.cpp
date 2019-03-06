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
#define WIDTH 5
#define HEIGHT 5
#define DEPTH 5
#define STEP 1

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

typedef Vec3<unsigned int> Vec3i;

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
    float transparency;
    float refractivity;
    Vec3f emissionColor;
};

struct light
{
    Vec3f direction;
    Vec3f color;
    float intensity;
};

float slope2 (float p1, float p2){
//    std::cout << p1 << ", " << p2 << std::endl;
//    std::cout << p2 - p1 << std::endl;
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
void inSphere(object *field,  Vec3f center, float radius, Vec3f surfaceColor, float reflectivity, float transparency, float refractivity, Vec3f emissionColor){
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
                field[i + j * WIDTH + k * WIDTH * HEIGHT].transparency = transparency;
                field[i + j * WIDTH + k * WIDTH * HEIGHT].emissionColor = emissionColor;
                field[i + j * WIDTH + k * WIDTH * HEIGHT].refractivity = refractivity;
            }
        }
    }
}

void calcLightSource(object *field, light light[], unsigned int direction, Vec3f color, float intensity){
    unsigned int coordinate;
    for (unsigned int i=0; i<WIDTH; ++i){
       for (unsigned int j=0; j<HEIGHT; ++j){
           for (unsigned int k=0; k<DEPTH; ++k){
               coordinate = i * (int)pow(WIDTH, direction) + j * (int)pow(WIDTH, (direction==0)?2:(direction - 1)) + k * (int)pow(WIDTH, (direction<2)?(direction + 1):(0));
               Vec3f previousColor;
               unsigned int previousCoo;
               if (i == 0){
                    previousColor = color;
               } else {
                   previousCoo = coordinate - (int)pow(WIDTH, direction);
                   previousColor = field[previousCoo].color;
               }
//               std::cout << i << ", " << j << ", " << k << ", " << coordinate << ", " << std::endl;
               light[coordinate].color = previousColor * field[coordinate].color * (1 - field[coordinate].transparency) * intensity;
           }
       }
    }
}

//[comment]
// This variable controls the maximum recursion depth
//[/comment]
#define MAX_RAY_DEPTH 5

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
Vec3f trace(
        const Vec3f &rayorig,
        const Vec3f &raydir,
        const std::vector<Sphere> &spheres,
        const int &depth)
{
    //if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
    float tnear = INFINITY;
    const Sphere* sphere = NULL;
    // find intersection of this ray with the sphere in the scene
    for (unsigned i = 0; i < spheres.size(); ++i) {
        float t0 = INFINITY, t1 = INFINITY;
        if (spheres[i].intersect(rayorig, raydir, t0, t1)) {
            if (t0 < 0) t0 = t1;
            if (t0 < tnear) {
                tnear = t0;
                sphere = &spheres[i];
            }
        }
    }
    // if there's no intersection return black or background color
    if (!sphere) return Vec3f(2);
    Vec3f surfaceColor = 0; // color of the ray/surfaceof the object intersected by the ray
    Vec3f phit = rayorig + raydir * tnear; // point of intersection
    Vec3f nhit = phit - sphere->center; // normal at the intersection point
    nhit.normalize(); // normalize normal direction
    // If the normal and the view direction are not opposite to each other
    // reverse the normal direction. That also means we are inside the sphere so set
    // the inside bool to true. Finally reverse the sign of IdotN which we want
    // positive.
    float bias = 1e-4; // add some bias to the point from which we will be tracing
    bool inside = false;
    if (raydir.dot(nhit) > 0) nhit = -nhit, inside = true;
    if ((sphere->transparency > 0 || sphere->reflection > 0) && depth < MAX_RAY_DEPTH) {
        float facingratio = -raydir.dot(nhit);
        // change the mix value to tweak the effect
        float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
        // compute reflection direction (not need to normalize because all vectors
        // are already normalized)
        Vec3f refldir = raydir - nhit * 2 * raydir.dot(nhit);
        refldir.normalize();
        Vec3f reflection = trace(phit + nhit * bias, refldir, spheres, depth + 1);
        Vec3f refraction = 0;
        // if the sphere is also transparent compute refraction ray (transmission)
        if (sphere->transparency) {
            float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
            float cosi = -nhit.dot(raydir);
            float k = 1 - eta * eta * (1 - cosi * cosi);
            Vec3f refrdir = raydir * eta + nhit * (eta *  cosi - sqrt(k));
            refrdir.normalize();
            refraction = trace(phit - nhit * bias, refrdir, spheres, depth + 1);
        }
        // the result is a mix of reflection and refraction (if the sphere is transparent)
        surfaceColor = (
                               reflection * fresneleffect +
                               refraction * (1 - fresneleffect) * sphere->transparency) * sphere->surfaceColor;
    }
    else {
        // it's a diffuse object, no need to raytrace any further
        for (unsigned i = 0; i < spheres.size(); ++i) {
            if (spheres[i].emissionColor.x > 0) {
                // this is a light
                Vec3f transmission = 1;
                Vec3f lightDirection = spheres[i].center - phit;
                lightDirection.normalize();
                for (unsigned j = 0; j < spheres.size(); ++j) {
                    if (i != j) {
                        float t0, t1;
                        if (spheres[j].intersect(phit + nhit * bias, lightDirection, t0, t1)) {
                            transmission = 0;
                            break;
                        }
                    }
                }
                surfaceColor += sphere->surfaceColor * transmission * std::max(float(0), nhit.dot(lightDirection)) * spheres[i].emissionColor;
            }
        }
    }

    return surfaceColor + sphere->emissionColor;
}

//[comment]
// Main rendering function. We compute a camera ray for each pixel of the image
// trace it and return a color. If the ray hits a sphere, we return the color of the
// sphere at the intersection point, else we return the background color.
//[/comment]
Vec3f *render(const std::vector<Sphere> &spheres)
{
    Vec3f *image = new Vec3f[WIDTH * HEIGHT], *pixel = image;
    float invWidth = 1 / float(WIDTH), invHeight = 1 / float(HEIGHT);
    float fov = 30, aspectratio = WIDTH / float(HEIGHT);
    float angle = tan(M_PI * 0.5 * fov / 180.);
    // Trace rays
    for (unsigned y = 0; y < HEIGHT; ++y) {
        for (unsigned x = 0; x < WIDTH; ++x, ++pixel) {
            float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
            float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
            Vec3f raydir(xx, yy, -1);
            raydir.normalize();
            *pixel = trace(Vec3f(0, 0, 0), raydir, spheres, 0);
        }
    }
    // Save result to a PPM image (keep these flags if you compile under Windows)
    std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);
    ofs << "P6\n" << WIDTH << " " << HEIGHT << "\n255\n";
    for (unsigned i = 0; i < WIDTH * HEIGHT; ++i) {
        ofs << (unsigned char) (std::min(float(1), image[i].x) * 255) <<
            (unsigned char) (std::min(float(1), image[i].y) * 255) <<
            (unsigned char) (std::min(float(1), image[i].z) * 255);
    }
    ofs.close();
//    delete [] image;
    return image;
//    return rgbImage;
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

    std::vector<Sphere> spheres;
    // 可変長配列vector、一つずつそれぞれがsphereクラス、push_backは末尾に追加
    // position, radius, surface color, reflectivity, transparency, emission color
    spheres.push_back(Sphere(Vec3f( 0.0, -10004, -20), 10000, Vec3f(0.20, 0.20, 0.20), 0, 0.0));
    spheres.push_back(Sphere(Vec3f( 0.0,      0, -20),     4, Vec3f(1.00, 0.32, 0.36), 1, 0.5));
    spheres.push_back(Sphere(Vec3f( 5.0,     -1, -15),     2, Vec3f(0.90, 0.76, 0.46), 1, 0.0));
    spheres.push_back(Sphere(Vec3f( 5.0,      0, -25),     3, Vec3f(0.65, 0.77, 0.97), 1, 0.0));
    spheres.push_back(Sphere(Vec3f(-5.5,      0, -15),     3, Vec3f(0.90, 0.90, 0.90), 1, 0.0));
    // light
    spheres.push_back(Sphere(Vec3f( -0.0,     200, -30),     3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(3)));

    std::cout << typeid(render(spheres)[2000].x).name() << std::endl;
    glDrawPixels(WIDTH, HEIGHT, GL_RGB, GL_FLOAT, render(spheres));
    glFlush();

//    glutSwapBuffers();
//    glDisable(GL_DEPTH_TEST);
//    glutSwapBuffers();
}

//void myReshape(int WIDTH, int HEIGHT)
//{
//    float aspect = (float) WIDTH / (float) HEIGHT;
//    glViewport(0, 0, WIDTH, HEIGHT);
//    glMatrixMode(GL_PROJECTION);
//    glLoadIdentity();
//    gluPerspective(30.0, aspect, 1.0, 50.0);
//    glMatrixMode(GL_MODELVIEW);
//}

//[comment]
// In the main function, we will create the scene which is composed of 5 spheres
// and 1 light (which is also a sphere). Then, once the scene description is complete
// we render that scene, by calling the render() function.
//[/comment]

int main(int argc, char **argv)
{
//    glutInit(&argc, argv);
//    myInit(argv[0]);
//    glutDisplayFunc(display);
//    glutReshapeFunc (myReshape);
//    glutMainLoop();

    std::cout << "start" << std::endl;
    object field[WIDTH * HEIGHT * DEPTH];
    light fromLightSource[WIDTH * HEIGHT * DEPTH];
    Vec3f grad[WIDTH * HEIGHT * DEPTH];
    for (unsigned int i=0; i < WIDTH * HEIGHT * DEPTH; ++i){
        field[i] = {0};
    }
//    field, center, radius, color, reflection, transparency, refraction, emissionColor
    inSphere(field,  Vec3f (0., 0., 0.), 2, Vec3f (1., 0., 0.), 0., 0., 0.5, Vec3f(0., 0., 0.));
    gradRefraction(field, grad);
//    for (unsigned int i=0; i<WIDTH * HEIGHT * DEPTH; ++i){
//        std::cout << grad[i] << std::endl;
//    }
//  field, light, direction, color, intensity
    calcLightSource(field, fromLightSource, 1, (1,1,1), 1.0);
//    std::cout << field[1].color << std::endl;
    return 0;
}
