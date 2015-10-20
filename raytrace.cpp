//
// template-rt.cpp
//

#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
using namespace std;

int g_width;
int g_height;

struct Ray
{
    vec4 origin;
    vec4 dir;
};

// TODO: add structs for spheres, lights and anything else you may need.

struct Sphere
{
    string name;
    vec3 pos;
    vec3 scale;
    vec4 color;
    vec4 K;
    float n;
};

struct Light
{
    string name;
    vec3 pos;
    vec3 I;
};

vector<vec4> g_colors;
vector<Sphere> spheres;
vector<Light> lights;
vec4 bg_color;
vec3 ambient;
string name;

float g_left;
float g_right;
float g_top;
float g_bottom;
float g_near;

vec4 getDir(int ix, int iy);
vec4 trace(const Ray& ray, int depth);

// -------------------------------------------------------------------
// Input file parsing

vec4 toVec4(const string& s1, const string& s2, const string& s3)
{
    stringstream ss(s1 + " " + s2 + " " + s3);
    vec4 result;
    ss >> result.x >> result.y >> result.z;
    result.w = 1.0f;
    return result;
}

float toFloat(const string& s)
{
    stringstream ss(s);
    float f;
    ss >> f;
    return f;
}

void parseLine(const vector<string>& vs)
{
    //TODO: add parsing of NEAR, LEFT, RIGHT, BOTTOM, TOP, SPHERE, LIGHT, BACK, AMBIENT, OUTPUT.
    if (vs[0] == "RES")
    {
        g_width = (int)toFloat(vs[1]);
        g_height = (int)toFloat(vs[2]);
        g_colors.resize(g_width * g_height);
    }
    
    if (vs[0] == "NEAR")
    {
        g_near = toFloat(vs[1]);
    }
    
    if (vs[0] == "LEFT")
    {
        g_left = toFloat(vs[1]);
    }
    
    if (vs[0] == "RIGHT")
    {
        g_right = toFloat(vs[1]);
    }
    
    if (vs[0] == "BOTTOM")
    {
        g_bottom = toFloat(vs[1]);
    }
    
    if (vs[0] == "TOP")
    {
        g_top = toFloat(vs[1]);
    }
    
    if (vs[0] == "SPHERE")
    {
        Sphere s;
        s.name = vs[1];
        s.pos = vec3(toFloat(vs[2]), toFloat(vs[3]), toFloat(vs[4]));
        s.scale = vec3(toFloat(vs[5]), toFloat(vs[6]), toFloat(vs[7]));
        s.color = vec4(toFloat(vs[8]), toFloat(vs[9]), toFloat(vs[10]), 1.0f);
        s.K = vec4(toFloat(vs[11]), toFloat(vs[12]), toFloat(vs[13]), toFloat(vs[14]));
        s.n = toFloat(vs[15]);
        spheres.push_back(s);
    }
    
    if (vs[0] == "LIGHT")
    {
        Light li;
        li.name = vs[1];
        li.pos = vec3(toFloat(vs[2]), toFloat(vs[3]), toFloat(vs[4]));
        li.I = vec3(toFloat(vs[5]), toFloat(vs[6]), toFloat(vs[7]));
        lights.push_back(li);
    }
    
    if (vs[0] == "BACK")
    {
        bg_color = vec4(toFloat(vs[1]), toFloat(vs[2]), toFloat(vs[3]), 0.0f);
    }
    
    if (vs[0] == "AMBIENT")
    {
        ambient = vec3(toFloat(vs[1]), toFloat(vs[2]), toFloat(vs[3]));
    }
    
    if (vs[0] == "OUTPUT")
    {
        name = vs[1];
    }
}

void loadFile(const char* filename)
{
    ifstream is(filename);
    if (is.fail())
    {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    string s;
    vector<string> vs;
    while(!is.eof())
    {
        vs.clear();
        getline(is, s);
        istringstream iss(s);
        while (!iss.eof())
        {
            string sub;
            iss >> sub;
            vs.push_back(sub);
        }
        parseLine(vs);
    }
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}


// -------------------------------------------------------------------
// Intersection routine

// TODO: add your ray-sphere intersection routine here.
bool findIntersection(const Ray& ray, float length)
{
    vec4 sval = ray.origin;
    vec4 cval = ray.dir;
    vec4 sprime;
    vec4 cprime;
    mat4 sphereMat;
    mat4 matInverse_forSphere;
    
    for (int i=0; i<spheres.size(); i++)
    {
        sphereMat=mat4(vec4(spheres[i].scale.x, 0.0f, 0.0f, spheres[i].pos.x),          // transformation matrix
                       vec4(0.0f, spheres[i].scale.y, 0.0f, spheres[i].pos.y),
                       vec4(0.0f, 0.0f, spheres[i].scale.z, spheres[i].pos.z),
                       vec4(0.0f, 0.0f, 0.0f, 1.0f));
        
        mat4 matInverse;                                                                // inverse matrix
        bool result = InvertMatrix(sphereMat, matInverse);
        if (result == false)
            printf("Matrix not invertible.\n");
        
        sprime = matInverse * sval;
        cprime = matInverse * cval;
        
        // quadratic values
        float clen = 1.*sqrt(1.*pow(cprime.x, 2) + 1.*pow(cprime.y, 2) + 1.*pow(cprime.z, 2));      // length of c
        float slen = 1.*sqrt(1.*pow(sprime.x, 2) + 1.*pow(sprime.y, 2) + 1.*pow(sprime.z, 2));      // length of s
        float product = sprime.x*cprime.x + sprime.y*cprime.y + sprime.z*cprime.z;                  // sval * cval
        float val = pow(product, 2) - (pow(clen, 2) * (pow(slen, 2) -1));                           // B^2 - AC
        float t1, t2;
        float tval;
        
        // find t
        if (val > 0)                                                                    // if there are two solutions for t
        {
            t1 = (-(1.*product)/(1.*pow(clen, 2))) + (1.*sqrt(val) / (1.*pow(clen, 2)));
            
            if (t1 > 1)
            {
                vec4 t1_intersect = ray.dir*t1;
                float t1_intersect_len = 1.*sqrt(1.*pow(t1_intersect.x,2)+1.*pow(t1_intersect.y,2)+1.*pow(t1_intersect.z,2));
            
                if(0.001f < t1 && t1 < length)
                {
                    return true;
                }
            }
            
            t2 = (-(1.*product)/(1.*pow(clen, 2))) - (1.*sqrt(val) / (1.*pow(clen, 2)));
            
            if (t2 > 1)
            {
                vec4 t2_intersect = ray.dir*t2;
                float t2_intersect_len = 1.*sqrt(1.*pow(t2_intersect.x,2)+1.*pow(t2_intersect.y,2)+1.*pow(t2_intersect.z,2));
            
                if(0.001f < t2 && t2 < length)
                {
                    return true;
                }
            }
        }
        if (val == 0)
        {
            tval = (1.*-product) / (1.*pow(clen, 2));
            
            if (tval > 1)
            {
                vec4 tval_intersect = ray.dir*tval;
                float tval_intersect_len = 1.*sqrt(1.*pow(tval_intersect.x,2)+1.*pow(tval_intersect.y,2)+1.*pow(tval_intersect.z,2));
                
                if (0.001f < tval && tval < length)
                {
                    return true;
                }
            }
        }
    }
    return false;
}

vec4 setPixelColor(int spherenum, vec4 sval, vec4 cval, vec4 sprime, vec4 cprime, float tmin, mat4 matInverse_forSphere, int depth)
{
    vec4 pixel_color(0.0f,0.0f,0.0f,1.0f);
    float rval = 0.0f;
    float gval = 0.0f;
    float bval = 0.0f;
    pixel_color = vec4(spheres[spherenum].K.x * ambient.x * spheres[spherenum].color.x,
                   spheres[spherenum].K.x * ambient.y * spheres[spherenum].color.y,
                   spheres[spherenum].K.x * ambient.z * spheres[spherenum].color.z,
                   1.0f);

    // Step 3: add illumination
    
    // hit point of transformed sphere
    vec4 P_trans = vec4(sval.x + 1.*(tmin * cval.x), sval.y + 1.*(tmin * cval.y), sval.z + 1.*(tmin * cval.z), 1.0f);

    // hit point of unit sphere
    vec4 P_unit;                                                                        // hit point of unit sphere
    P_unit = vec4(sprime.x + 1.*(tmin * cprime.x), sprime.y + 1.*(tmin * cprime.y), sprime.z + 1.*(tmin * cprime.z), 1.0f);

    // unit normal vector
    vec4 unit_normal = vec4(sprime.x + 1.*(tmin * cprime.x), sprime.y + 1.*(tmin * cprime.y), sprime.z + 1.*(tmin * cprime.z), 0.0f);

    // transformed normal vector
    vec4 trans_normal = transpose(matInverse_forSphere) * unit_normal;
    float trans_normallen=1.*sqrt(1.*pow(1.*trans_normal.x,2)+1.*pow(1.*trans_normal.y,2)+1.*pow(1.*trans_normal.z,2));
    vec4 new_trans=vec4((1.*trans_normal.x)/(1.*trans_normallen),(1.*trans_normal.y)/(1.*trans_normallen),(1.*trans_normal.z)/(1.*trans_normallen),0.0f);

    // reflection
    Ray reflect;
    reflect.origin = vec4(sval.x + 1.*(tmin * cval.x), sval.y + 1.*(tmin * cval.y), sval.z + 1.*(tmin * cval.z), 1.0f);
    
    // v = -2(N dot c)N + c
    vec4 v = -2 * 1.*dot(new_trans, cval) * 1.*new_trans + cval;
    float v_len = 1.*sqrt(1.*pow(v.x,2)+1.*pow(v.y,2)+1.*pow(v.z,2));
    vec4 new_v = vec4(1.*v.x/(1.*v_len), 1.*v.y/(1.*v_len), 1.*v.z/(1.*v_len), 0.0f);
    reflect.dir = new_v;
    vec4 reflect_color = trace(reflect, depth);
    if(reflect_color.x == bg_color.x && reflect_color.y == bg_color.y && reflect_color.z == bg_color.z)
        reflect_color = vec4(0.0f,0.0f,0.0f,1.0f);
    reflect_color *= spheres[spherenum].K.w;
    
    bool shadow=false;
    for (int j=0; j<lights.size(); j++)                                                 // for each light source
    {
        // L
        vec4 L = vec4(1.*(lights[j].pos.x - P_trans.x), 1.*(lights[j].pos.y - P_trans.y), 1.*(lights[j].pos.z - P_trans.z), 0.0f);
        float L_len = 1.*sqrt(1.*pow(1.*L.x, 2) + 1.*pow(1.*L.y, 2) + 1.*pow(1.*L.z, 2));
        vec4 newL = vec4((1.*L.x)/(1.*L_len), (1.*L.y)/(1.*L_len), (1.*L.z)/(1.*L_len), 0.0f);
        
        // r
        vec4 reflect_vec = (1.*(2*new_trans)*(1.*dot(new_trans, newL)) - newL);
        float reflect_veclen = 1.*sqrt(1.*pow(reflect_vec.x, 2)+1.*pow(reflect_vec.y, 2)+1.*pow(reflect_vec.z, 2));
        vec4 new_reflect = vec4(1.*reflect_vec.x/(1.*reflect_veclen), 1.*reflect_vec.y/(1.*reflect_veclen), 1.*reflect_vec.z/(1.*reflect_veclen), 0.0f);
        
        // v
        vec4 eye_vec = sval - P_trans;
        float eye_veclen = 1.*sqrt(1.*pow(eye_vec.x, 2)+1.*pow(eye_vec.y, 2)+1.*pow(eye_vec.z, 2));
        vec4 new_eye = vec4(1.*eye_vec.x/(1.*eye_veclen), 1.*eye_vec.y/(1.*eye_veclen), 1.*eye_vec.z/(1.*eye_veclen), 0.0f);
        
        // account for shadowing
        Ray ray;
        
        vec4 c = vec4(1.*lights[j].pos.x - 1.*P_trans.x, 1.*lights[j].pos.y - 1.*P_trans.y, 1.*lights[j].pos.z - 1.*P_trans.z, 0.0f);
        float c_len = sqrt(1.*pow(c.x,2)+1.*pow(c.y,2)+1.*pow(c.z,2));
        vec4 new_c = vec4(1.*c.x/(1.*c_len), 1.*c.y/(1.*c_len), 1.*c.z/(1.*c_len), 0.0f);
        
        ray.origin = P_trans;
        ray.dir = new_c;
        
        shadow = findIntersection(ray, c_len);
        
        if (shadow == false)
        {
            if(dot(new_trans,newL)>=0)
            {
                // add diffuse components
                float K_diff = 1.*spheres[spherenum].K.y;
                rval += (1.*lights[j].I.x * 1.*K_diff * (1.*dot(new_trans, newL)) * spheres[spherenum].color.x);
                gval += (1.*lights[j].I.y * 1.*K_diff * (1.*dot(new_trans, newL)) * spheres[spherenum].color.y);
                bval += (1.*lights[j].I.z * 1.*K_diff * (1.*dot(new_trans, newL)) * spheres[spherenum].color.z);
            
                if (dot(new_reflect, new_eye) >= 0)
                {
                    // add specular components
                    float K_spec = 1.*spheres[spherenum].K.z;
                    rval += K_spec * 1.*lights[j].I.x * (1.*pow((1.*dot(new_reflect, new_eye)), spheres[spherenum].n));
                    gval += K_spec * 1.*lights[j].I.y * (1.*pow((1.*dot(new_reflect, new_eye)), spheres[spherenum].n));
                    bval += K_spec * 1.*lights[j].I.z * (1.*pow((1.*dot(new_reflect, new_eye)), spheres[spherenum].n));
                }
            }
        }
    }
    
    

    pixel_color.x += rval;
    pixel_color.y += gval;
    pixel_color.z += bval;

    
    
    return pixel_color + reflect_color;
}

// -------------------------------------------------------------------
// Ray tracing

vec4 trace(const Ray& ray, int depth)
{
    // TODO: implement your ray tracing routine here.
    
    // Step 2: no illumination
    vec4 sval = ray.origin;
    vec4 cval = ray.dir;
    
    float tmin;             // smallest t
    bool hasTval = true;
    bool firstT = true;
    int spherenum = 0;      // index of the sphere with the smallest t
    int noTval = 0;
    vec4 sprime;
    vec4 cprime;
    mat4 sphereMat;
    mat4 matInverse_forSphere;
    vec4 pixel_color;
    
    
    if(depth==3)
        return pixel_color;
    depth++;
    for (int i=0; i<spheres.size(); i++)
    {
        sphereMat=mat4(vec4(spheres[i].scale.x, 0.0f, 0.0f, spheres[i].pos.x),          // transformation matrix
                       vec4(0.0f, spheres[i].scale.y, 0.0f, spheres[i].pos.y),
                       vec4(0.0f, 0.0f, spheres[i].scale.z, spheres[i].pos.z),
                       vec4(0.0f, 0.0f, 0.0f, 1.0f));
        
        mat4 matInverse;                                                                // inverse matrix
        bool result = InvertMatrix(sphereMat, matInverse);
        if (result == false)
            printf("Matrix not invertible.\n");
        
        sprime = matInverse * sval;
        cprime = matInverse * cval;
        
        float clen = 1.*sqrt(1.*pow(cprime.x, 2) + 1.*pow(cprime.y, 2) + 1.*pow(cprime.z, 2));      // length of c
        float slen = 1.*sqrt(1.*pow(sprime.x, 2) + 1.*pow(sprime.y, 2) + 1.*pow(sprime.z, 2));      // length of s
        float product = sprime.x*cprime.x + sprime.y*cprime.y + sprime.z*cprime.z;                  // sval * cval
        float val = pow(product, 2) - (pow(clen, 2) * (pow(slen, 2) -1));                           // B^2 - AC
        float t1, t2;
        float tval;
       
        // find t
        if (val > 0)                                                                    // if there are two solutions for t
        {
            t1 = (-(1.*product)/(1.*pow(clen, 2))) + (1.*sqrt(val) / (1.*pow(clen, 2)));
            t2 = (-(1.*product)/(1.*pow(clen, 2))) - (1.*sqrt(val) / (1.*pow(clen, 2)));
        
            if (t1 > 1||(t1>0.0001&&depth>1))
            {
                if(firstT)                                                              // set only the first t to tmin
                {
                    tmin = t1;
                    firstT = false;
                    spherenum = i;
                    matInverse_forSphere = matInverse;
                    pixel_color=setPixelColor(spherenum, sval, cval, sprime, cprime, tmin, matInverse_forSphere, depth);
                }
                if(t1 < tmin)                                                           // compare t to tmin
                {
                    tmin = t1;
                    spherenum = i;
                    matInverse_forSphere = matInverse;
                    pixel_color=setPixelColor(spherenum, sval, cval, sprime, cprime, tmin, matInverse_forSphere, depth);
                }
            }
            if (t2 > 1||(t2>0.0001&&depth>1))
            {
                if(firstT)
                {
                    tmin = t2;
                    firstT = false;
                    spherenum = i;
                    matInverse_forSphere = matInverse;
                    pixel_color=setPixelColor(spherenum, sval, cval, sprime, cprime, tmin, matInverse_forSphere, depth);
                }
                if (t2 < tmin)
                {
                    tmin = t2;
                    spherenum = i;
                    matInverse_forSphere = matInverse;
                    pixel_color=setPixelColor(spherenum, sval, cval, sprime, cprime, tmin, matInverse_forSphere, depth);
                }
            }
            else if (t1 <= 1 && t2 <= 1)
            {
                noTval++;
            }
        }
        
        else if (val == 0)                                                              // if there is one solution for t
        {
            tval = (1.*-product) / (1.*pow(clen, 2));
            if (tval > 1||(tval>0.0001&&depth>1))
            {
                if(firstT)
                {
                    tmin = tval;
                    firstT = false;
                    spherenum = i;
                    matInverse_forSphere = matInverse;
                    pixel_color=setPixelColor(spherenum, sval, cval, sprime, cprime, tmin, matInverse_forSphere, depth);
                }
                if (tval < tmin)
                {
                    tmin = tval;
                    spherenum = i;
                    matInverse_forSphere = matInverse;
                    pixel_color=setPixelColor(spherenum, sval, cval, sprime, cprime, tmin, matInverse_forSphere, depth);
                }
            }
        }
        
        else if (val < 0)                                                               // no solutions for t
        {
            noTval++;
        }
        
        if (noTval >= spheres.size())                                                   // number of intersection fails > number of spheres
        {                                                                               // means there is no t value
            hasTval = false;
        }
    }
    
    if (!hasTval)                                                                       // no t value, set color to background color
    {
        pixel_color = vec4(bg_color.x, bg_color.y, bg_color.z, 1.0f);
    }
    
    // reflection
   
    
    return pixel_color;
    //return vec4(0.0f, 0.0f, 0.0f, 1.0f);
}

vec4 getDir(int ix, int iy)
{
    // TODO: modify this. This should return the direction from the origin
    // to pixel (ix, iy), normalized.
    vec4 dir;
    float x = g_left + (1.*ix/(1.*g_width))*(g_right - g_left);
    float y = g_bottom + (1.*iy/(1.*g_height))*(g_top - g_bottom);
    float z = -g_near;
    dir = vec4(x, y, z, 0.0f);
    //dir = vec4(0.0f, 0.0f, -1.0f, 0.0f);
    return dir;
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    ray.dir = getDir(ix, iy);
    int depth = 0;
    vec4 color = trace(ray, depth);
    setColor(ix, iy, color);
}

void render()
{
    for (int iy = 0; iy < g_height; iy++)
        for (int ix = 0; ix < g_width; ix++)
            renderPixel(ix, iy);
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, const char* fname, unsigned char* pixels)
{
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }

    fclose(fp);
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    // TODO: clamp values if out of range.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++)
            {
                if (g_colors[y*g_width + x][i] > 1.0) { g_colors[y*g_width + x][i] = 1.0f; }
                if (g_colors[y*g_width + x][i] < 0.0) { g_colors[y*g_width + x][i] = 0.0f; }
                buf[y*g_width*3+x*3+i] = (unsigned char)(((float*)g_colors[y*g_width+x])[i] * 255.9f);
            }
    
    // TODO: change file name based on input file name.
    savePPM(g_width, g_height, name.c_str(), buf);
    delete[] buf;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: template-rt <input_file.txt>" << endl;
        exit(1);
    }
    loadFile(argv[1]);
    render();
    saveFile();
	return 0;
}

