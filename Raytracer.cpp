//Brian Pingel
//COEN 148
//Raytracing and Shaders
/*I have gotten most of the calculations correct
Part 1: Phong Shading works exactly as intended
Part 2: Reflections works perfectly as well
Part 3: Refractions working partially,
The sphere is transparent and has a reflection in it but the refraction does not work properly so refraction is not entirely functional, there is refraction like behavior however there is a math error somewhere in the
    calculation of the refraction vector*/


#include <stdio.h>
#ifndef _WIN32
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <math.h>
#include <iostream>


using namespace std;


//Improves efficiency of dividing by a square root which is inefficient due to floating point division being slow
//Taken from Quake III Arena Source Code
float Q_rsqrt(float number)
{
        long i;
        float x2, y;
        const float threehalfs = 1.5F;


        x2 = number * 0.5F;
        y = number;
        i = *(long *)&y;                       // evil floating point bit level hacking
        i = 0x5f3759df - (i >> 1); 
        y = *(float *)&i;
        y = y * (threehalfs - (x2 * y * y));   // 1st iteration
        //y = y * (threehalfs - (x2 * y * y));   // 2nd iteration, this can be removed


        return y;
}


class point {
public:
        float x;
        float y;
        float z;
        point(float initx = 0, float inity = 0, float initz = 0) {
                x = initx;
                y = inity;
                z = initz;
        }
        point& operator=(point& init) {
                x = init.x;
                y = init.y;
                z = init.z;
                return *this;
        }
        point& operator-=(point rhs) {
                float new_x = x - rhs.x;
                float new_y = y - rhs.y;
                float new_z = z - rhs.z;
                x = new_x;
                y = new_y;
                z = new_z;
                return *this;
        }
};


point operator-(point lhs, const point& rhs) {
        return lhs -= rhs;
}


class light {
public:
        point center;
        float intensity;
        float specular;
        light(float inx, float iny, float inz, float intensityIn, float specularIn) {
                intensity = intensityIn;
                specular = specularIn;
                center.x = inx;
                center.y = iny;
                center.z = inz;


        }
};


class vec {
public:
        vec(float initx = 0, float inity = 0, float initz = 0) : x(initx), y(inity), z(initz) {}
        vec(point data) : x(data.x), y(data.y), z(data.z) {}
        vec(point origin, point end) {
                x = (end.x - origin.x);
                y = (end.y - origin.y);
                z = (end.z - origin.z);
        }
        float x;
        float y;
        float z;
        void normalize() {
                float length = Q_rsqrt((x*x) + (y*y) + (z*z));
                x *= length;
                y *= length;
                z *= length;
        }
        friend double dot(vec, vec);
        friend vec cross(vec, vec);
        void mulvec(float mul) {
                x *= mul;
                y *= mul;
                z *= mul;
        }
        vec& operator=(vec rhs) {
                this->x = rhs.x;
                this->y = rhs.y;
                this->z = rhs.z;
                return *this;
        }
        vec& operator-=(vec rhs) {
                this->x -= rhs.x;
                this->y -= rhs.y;
                this->z -= rhs.z;
                return *this;
        }
};


double dot(vec vector1, vec vector2) {
        return ((vector1.x*vector2.x) + (vector1.y * vector2.y) + (vector1.z * vector2.z));
}
vec cross(vec a, vec b) {
        vec crossProd;
        crossProd.x = ((a.y * b.z) - (a.z * b.y));
        crossProd.y = ((a.z * b.x) - (a.x * b.z));
        crossProd.z = ((a.x * b.y) - (a.y * b.x));
        return crossProd;
}


class ray {
public:
        point origin;
        vec dir;
        ray(point begin, vec vector) {
                dir.x = vector.x;
                dir.y = vector.y;
                dir.z = vector.z;
                origin.x = begin.x;
                origin.y = begin.y;
                origin.z = begin.z;
        }
        ray(float x0, float y0, float z0, float xv, float yv, float zv) {
                origin.x = x0;
                origin.y = y0;
                origin.z = z0;
                dir.x = xv;
                dir.y = yv;
                dir.z = zv;
        }
};


class shape {
public:
        point center;
        double R;
        double G;
        double B;
        float Reflec;
        float Spec;
        int type;
        //initialize to the input xyz coords
        shape(float inxc = 0, float inyc = 0, float inzc = 0, double inR = 0, double inG = 0, double inB = 0, float new_reflec = 0, float new_illum = 0, int typeIn = 1) : R(inR), G(inG), B(inB), Reflec(new_reflec), Spec(new_illum), type(typeIn) {
                center.x = inxc;
                center.y = inyc;
                center.z = inzc;
        }
        virtual  bool intersect(point*, ray, double*) = 0;
        virtual vec normal(point) = 0;
};


class sphere : public shape {
public:
        float radius;
        sphere(float rad, float inxc = 0, float inyc = 0, float inzc = 0, double inR = 0, double inG = 0, double inB = 0, float new_reflec = 0, float new_illum = 0, int typeIn = 1) : shape(inxc, inyc, inzc, inR, inG, inB, new_reflec, new_illum, typeIn) {
                radius = rad;
        }
        //we utilize the discriminant to determine whether there is an intersection with a given sphere
        bool intersect(point *inter, ray tracer, double *t) {
                float rsq = radius * radius;
                tracer.dir.normalize();
                vec ominc(this->center, tracer.origin);
                double b = (dot(tracer.dir, ominc));
                double c = dot(ominc, ominc) - rsq;
                double discr = (b * b) - c;
                if (discr < 0) {
                        //no circ
                        return false;
                }
                else if (discr == 0) {
                        double d = -b;
                        tracer.dir.mulvec(d);
                        inter->x = tracer.dir.x + tracer.origin.x;
                        inter->y = tracer.dir.y + tracer.origin.y;
                        inter->z = tracer.dir.z + tracer.origin.z;
                        return true;
                        //exactly one intersect
                }


                else if (discr > 0) {
                        //multiple found so which is closer
                        double d1 = -b + sqrt(discr);
                        double d2 = -b - sqrt(discr);
                        if (d1 < d2 && d1 < *t && d1 > 0) {
                                *t = d1;
                                tracer.dir.mulvec(d1);
                                inter->x = tracer.dir.x + tracer.origin.x;
                                inter->y = tracer.dir.y + tracer.origin.y;
                                inter->z = tracer.dir.z + tracer.origin.z;
                        }
                        else if (d2 < d1 && d2 < *t && d2 > 0) {
                                *t = d2;
                                tracer.dir.mulvec(d2);
                                inter->x = tracer.dir.x + tracer.origin.x;
                                inter->y = tracer.dir.y + tracer.origin.y;
                                inter->z = tracer.dir.z + tracer.origin.z;
                        }
                        return true;
                }
                return false;
        }
        vec normal(point surface) {
                vec normalret(this->center, surface);
                return normalret;
        }
};


class plane : public shape {
public:
        vec upNormal;
        plane(double inx = 0, double iny = 0, double inz = 0, double R = .5, double G = .5, double B = .5, double reflecSet = .1, double SpecSet = .1, double typeIn = 1)
                : shape(inx, iny, inz, R, G, B, reflecSet, SpecSet, typeIn) {
                upNormal.x = 0;
                upNormal.z = 0;
                upNormal.y = 1;
        }
        bool intersect(point *inter, ray tracer, double *t) {
                point p;
                //calculates the distance variable
                p = center;
                p -= tracer.origin;
                double d = (dot(p, upNormal) / dot(upNormal, tracer.dir));
                if (d < 0)
                        return false;
                else if (d < *t) {
                        *t = d;
                        tracer.dir.mulvec(d);
                        //calculates the intersection point through the parameterization of the ray
                        inter->x = tracer.dir.x + tracer.origin.x;
                        inter->y = tracer.dir.y + tracer.origin.y;
                        inter->z = tracer.dir.z + tracer.origin.z;
                        return true;
                }
                return false;
        }
        vec normal(point surface) {
                return upNormal;
        }
};


class RGBColor {
public:
        double Red;
        double Green;
        double Blue;
        RGBColor(double Redin = 0, double Greenin = 0, double Bluein = 0) : Red(Redin), Green(Greenin), Blue(Bluein) {}
        RGBColor& operator=(const RGBColor rhs) {
                Red = rhs.Red;
                Green = rhs.Green;
                Blue = rhs.Blue;
                return *this;
        }
};


//calculates the max of two numbers
double max(double x1, double x2) {
        if (x1 > x2)
                return x1;
        else
                return x2;
}
const double windowHeighth = 1000;
const double windowWidth = 1000;
double windowSize = windowHeighth * windowWidth * 3;
point camera(3.5, 2, 0);
//light source
light source(3, 7, -4, .5, .9);
GLfloat pixels[1000 * 1000 * 3];
point viewPlaneCenter(2, 2, 0);
double planeWidth = 4;
double planeHeight = 4;
double planeX = 2;
double subdivs = planeWidth / windowWidth;
int maxDepth = 4;
double Airn = 1;
double Glassn = 1.52;


//objects in space
sphere tester(1, 0, 1, 0, .1, .8, .6, .8, .9, 1);
sphere tester2(1, 0, 3, -1, .7, .6, .9, .9, 1, 2);
sphere tester3(2, -2, 3, 4, .6, .9, .7, .6, 1, 3);
//sphere tester4(1.3, -2, 3, 4, .2, .1, .7, .4, 1, 2);
plane ground(0, 0, 0, .8, .5, .1, .9, .7, 1);
shape *objects[10];
int numObjects = 4;
double universalAmbient = .5;
double universalSpecular = 20;


/* raytrace receives a ray the recursive depth and a pointer to the pixel color to modify*/
void raytrace(ray tracer, int depth, RGBColor *pixelColor) {
        /*return case if gone too far, don't really care at a certain bounce
        more bounces means more time rendering a scene*/
        if (depth > maxDepth) {
                return;
        }
        //Part 1
        //variables to store object intersections
        int objectLoc = 0;
        bool foundOb = false;
        double storeT = 20000;
        double t = 20000000;
        point intersecter;
        point *interpoint;
        interpoint = &intersecter;
        int l = 0;
        //iterates through objects to determine intersections with ray and saves the intersected object
        for (l = 0; l < numObjects; l++) {
                foundOb = objects[l]->intersect(interpoint, tracer, &t);
                if (t < storeT && t > 0) {
                        storeT = t;
                        objectLoc = l;
                }
        }
        //if no object found
        if (!foundOb && l == numObjects) {
                pixelColor->Red *= .6;
                pixelColor->Green *= .6;
                pixelColor->Blue *= .6;
        }


        //calculate all necessary reflection and direction vectors of object for shading
        vec tolight(intersecter, source.center);
        tolight.normalize();
        vec normal;
        vec toCamera(intersecter, tracer.origin);
        toCamera.normalize();
        normal = objects[objectLoc]->normal(intersecter);
        normal.normalize();
        double twoDotLN = 2 * dot(tolight, normal);
        
        //calculate the perfect light reflection
        vec lightReflec = normal;
        lightReflec.mulvec(twoDotLN);
        lightReflec -= tolight;
        lightReflec.normalize();
        double costheta = dot(normal, tolight);
        //variables to calculate shadow intersects
        ray shadow(intersecter, tolight);
        point surrogate(0, 0, 0);
        point* shadowTest;
        shadowTest = &surrogate;
        bool shadowIntersec = false;
        //calculate diffuse, ambient, and specular contributions
        double Ispec = objects[objectLoc]->Spec * pow(max(0, dot(lightReflec, toCamera)), universalSpecular) * source.specular;
        double Idiff = objects[objectLoc]->Reflec * source.intensity *  max(0, costheta);
        double Iamb = objects[objectLoc]->Reflec * universalAmbient;
        //uses math to calculate the visual reflection
        vec eyeReflec(tracer.origin, intersecter);
        double dotprod = 2 * dot(eyeReflec, normal);
        vec normCopy = normal;
        normCopy.mulvec(dotprod);
        eyeReflec -= normCopy;
        eyeReflec.normalize();
        RGBColor recurseReflec(0, 0, 0);
        RGBColor recurseRefrac(0, 0, 0);
        RGBColor *RecurseReflecpoint = &recurseReflec;
        RGBColor *RecurseRefracpoint = &recurseRefrac;
        //determines whether or not there's a object blocking the path to the light source
        //if one is found then the color is darkened to look like a shadow
        for (int o = 0; o < numObjects; o++) {
                double p = 20902;
                shadowIntersec = objects[o]->intersect(shadowTest, shadow, &p);
                if (shadowIntersec) {
                        if (o == objectLoc)
                                continue;
                        else if (objects[objectLoc]->type == 1) {
                                //color reduction
                                pixelColor->Red = objects[objectLoc]->R / 10;
                                pixelColor->Green = objects[objectLoc]->G / 10;
                                pixelColor->Blue = objects[objectLoc]->B / 10;
                                return;
                        }
                        //Part 2
                        else if (objects[objectLoc]->type == 2 || objects[objectLoc]->type == 3) {
                                //even in the shadows a mirror can reflect;
                                ray reflectRay(intersecter, eyeReflec);
                                raytrace(reflectRay, depth + 1, RecurseReflecpoint);
                                //Part 3
                                double Rreflec = 1;
                                if (objects[objectLoc]->type == 3) {
                                        //determine refractive index ratios
                                        double nDiv = 0;
                                        double N1 = Airn;
                                        double N2 = Glassn;
                                        //Calculate refraction vector with Snell's law in vector form
                                        if (depth % 2 == 0) {
                                                nDiv = Airn / Glassn;
                                        }
                                        else {
                                                nDiv = Glassn / Airn;
                                                double ntemp = N1;
                                                N1 = N2;
                                                N2 = ntemp;
                                        }
                                        //depth = 0;
                                        double NormdotInc = dot(toCamera, normal);
                                        double nDivSq = nDiv * nDiv;
                                        double Sqroot = sqrt(1 - (nDivSq * (1 - NormdotInc)));
                                        vec NormalRefracCopy = normal;
                                        NormalRefracCopy.mulvec(Sqroot);
                                        vec NegNormal = normal;
                                        NegNormal.mulvec(-1);
                                        double m = nDiv * NormdotInc;
                                        vec crossNormS1 = cross(normal, toCamera);
                                        vec eyeRefrac = cross(NegNormal, crossNormS1);
                                        eyeRefrac.mulvec(nDiv);
                                        eyeRefrac -= NormalRefracCopy;
                                        ray refracRay(intersecter, eyeRefrac);
                                        raytrace(refracRay, depth + 1, RecurseRefracpoint);
                                        vec InNormal = normal;
                                        InNormal.mulvec(-1);
                                        //calculate Fresnel Effect Values
                                        double RefracAngle = dot(InNormal, eyeRefrac);
                                        double Rs = ((N1 * NormdotInc) - (N2 * RefracAngle))/((N1 * NormdotInc) + (N2 * RefracAngle));
                                        Rs *= Rs;
                                        double Rp = ((N2 * RefracAngle) - (N1 * NormdotInc)) / ((N2 * RefracAngle) + (N1 * NormdotInc));
                                        Rp *= Rp;
                                        Rreflec = (Rs + Rp) / 2;
                                }
                                double Rrefrac = 1 - Rreflec;
                                pixelColor->Red = objects[objectLoc]->Reflec * ((objects[objectLoc]->R) * (Ispec + (recurseReflec.Red * Rreflec) + (recurseRefrac.Red * Rrefrac)));
                                pixelColor->Green = objects[objectLoc]->Reflec * ((objects[objectLoc]->G) * (Ispec + (recurseReflec.Green * Rreflec) + (recurseRefrac.Green * Rrefrac)));
                                pixelColor->Blue = objects[objectLoc]->Reflec * ((objects[objectLoc]->B) * (Ispec + (recurseReflec.Blue * Rreflec) + (recurseRefrac.Blue * Rrefrac)));
                                return;
                        }
                }
        }
        if (objects[objectLoc]->type == 1) {
                //Specular lighting works, just a little hard to see with diffuse turned on as well
                //calculates phong shading
                pixelColor->Red = objects[objectLoc]->R * (Iamb + Idiff + Ispec);
                pixelColor->Green = objects[objectLoc]->G *  (Iamb + Idiff + Ispec);
                pixelColor->Blue = objects[objectLoc]->B * (Iamb + Idiff + Ispec);
                return;
        }
        else if (objects[objectLoc]->type == 2 || objects[objectLoc]->type == 3) {
                //Part 2
                //type 2 or 3 specifies at least reflective so it calls the recursive call for reflection
                ray reflectRay(intersecter, eyeReflec);
                raytrace(reflectRay, depth + 1, RecurseReflecpoint);
                double Rreflec = 1;
                if (objects[objectLoc]->type == 3) {
                        //Part 3
                        //type 3 specifies reflective refractive object
                        //determine refractive index ratios
                        double nDiv = 0;
                        double N1 = Airn;
                        double N2 = Glassn;
                        if (depth % 2 == 0) {
                                nDiv = Airn / Glassn;
                        }
                        else {
                                nDiv = Glassn / Airn;
                                double ntemp = N1;
                                N1 = N2;
                                N2 = ntemp;
                        }
                        //depth = 0;
                        //calculate the refraction vector withh Snell's Law in vector form;
                        vec toCameraCopy = toCamera;
                        double NormdotInc = dot(toCameraCopy, normal);
                        double nDivSq = nDiv * nDiv;
                        double Sqroot = sqrt(1 - (nDivSq * (1 - NormdotInc)));
                        vec NormalRefracCopy = normal;
                        NormalRefracCopy.mulvec(Sqroot);
                        vec NegNormal = normal;
                        NegNormal.mulvec(-1);
                        double m = nDiv * NormdotInc;
                        vec crossNormS1 = cross(NegNormal, toCamera);
                        vec eyeRefrac = cross(normal, crossNormS1);
                        eyeRefrac.mulvec(nDiv);
                        eyeRefrac -= NormalRefracCopy;
                        ray refracRay(intersecter, eyeRefrac);
                        raytrace(refracRay, depth + 1, RecurseRefracpoint);
                        vec InNormal = normal;
                        InNormal.mulvec(-1);
                        double RefracAngle = dot(InNormal, eyeRefrac);
                        double Rs = ((N1 * NormdotInc) - (N2 * RefracAngle)) / ((N1 * NormdotInc) + (N2 * RefracAngle));
                        Rs *= Rs;
                        double Rp = ((N2 * NormdotInc) - (N1 * dot(normal, eyeReflec))) / ((N2 * NormdotInc) + (N1 * RefracAngle));
                        Rp *= Rp;
                        Rreflec = (Rs + Rp) / 2;
                }
                double Rrefrac = 1 - Rreflec;
                pixelColor->Red = objects[objectLoc]->R * ((objects[objectLoc]->Reflec) * (Ispec + ((recurseReflec.Red * Rreflec) + (recurseRefrac.Red * Rrefrac))));
                pixelColor->Green = objects[objectLoc]->G * ((objects[objectLoc]->Reflec) * (Ispec + ((recurseReflec.Green * Rreflec) + (recurseRefrac.Green * Rrefrac))));
                pixelColor->Blue = objects[objectLoc]->B * ((objects[objectLoc]->Reflec) * (Ispec + ((recurseReflec.Blue * Rreflec) + (recurseRefrac.Green * Rrefrac))));
                return;
        }
}


void display(void) {
        //Print OpenGL errors 
        GLenum err_code;
        do {
                err_code = glGetError();
                if (err_code != GL_NO_ERROR)
                        printf("Error: %s\n", gluErrorString(err_code));
        } while (err_code != GL_NO_ERROR);
        objects[2] = &tester;
        objects[1] = &tester2;
        objects[0] = &ground;
        objects[3] = &tester3;
        //objects[4] = &tester4;


        //Clear buffer data 
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


        int i = 0;
        double j, z;
        //calculates vector from camera to specific point on the viewing plane
        for (j = 0; j <= 4; j += subdivs) {
                for (z = -2; z <= 2; i += 3, z += subdivs) {
                        point planeLoc(planeX, j, z);
                        vec rayvec(camera, planeLoc);
                        rayvec.normalize();
                        ray fire(camera, rayvec);
                        RGBColor* pixel;
                        RGBColor pixelColor(1, 1, 1);
                        pixel = &pixelColor;
                        raytrace(fire, 0, pixel);
                        //puts the pixelcolor variable into the pixels array
                        pixels[i] = pixelColor.Red;
                        pixels[i + 1] = pixelColor.Green;
                        pixels[i + 2] = pixelColor.Blue;
                }
        }
        cout << "done" << endl;
        glDrawPixels(windowHeighth, windowWidth, GL_RGB, GL_FLOAT, pixels);
        //Flush data 
        glFlush();
}


void init(void) {
        //Initialize lighting 
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_COLOR_MATERIAL);


        //Initialize camera 
        glMatrixMode(GL_PROJECTION);
        gluPerspective(50, 1, 0.1, 10);
        glMatrixMode(GL_MODELVIEW);


}


int main(int argc, char** argv) {
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);


        glutInitWindowSize(windowHeighth, windowWidth);
        glutInitWindowPosition(100, 100);
        glutCreateWindow("Assignment 3: Coen 148: Brian Pingel");
        init();
        glutDisplayFunc(display);




        glutMainLoop();
        return 0;
}
