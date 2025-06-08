#pragma once
#include "body.h"
#include <vector>

#ifndef BH_H
#define BH_H

struct Quad
{
    double cx, cy; // center of quad
    double length;

    Quad(double cx, double cy, double length);
    bool contains(const Body &b) const;
    Quad NW() const;
    Quad NE() const;
    Quad SW() const;
    Quad SE() const;
};

class BHTree
{
public:
    Quad region;
    Body *body = nullptr;
    double mass = 0.0;
    double com_x = 0.0;
    double com_y = 0.0;
    BHTree *NW = nullptr;
    BHTree *NE = nullptr;
    BHTree *SW = nullptr;
    BHTree *SE = nullptr;

    BHTree(Quad region);
    ~BHTree();

    void insert(Body *b);
    void update_force(Body *b, double theta, double G);
};

#endif // BH_H