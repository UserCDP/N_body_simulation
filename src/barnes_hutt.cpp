#include "barnes_hutt.h"
#include <cmath>

Quad::Quad(double cx, double cy, double length) : cx(cx), cy(cy), length(length) {}
bool Quad::contains(const Body &b) const
{
    return std::abs(b.position[0] - cx) <= length / 2 &&
           std::abs(b.position[1] - cy) <= length / 2;
}
Quad Quad::NW() const { return Quad(cx - length / 4, cy + length / 4, length / 2); }
Quad Quad::NE() const { return Quad(cx + length / 4, cy + length / 4, length / 2); }
Quad Quad::SW() const { return Quad(cx - length / 4, cy - length / 4, length / 2); }
Quad Quad::SE() const { return Quad(cx + length / 4, cy - length / 4, length / 2); }

// Barnes-Hut Tree Implementation
BHTree::BHTree(Quad region) : region(region) {}
BHTree::~BHTree()
{
    delete NW;
    delete NE;
    delete SW;
    delete SE;
}

void BHTree::insert(Body *b)
{
    if (!region.contains(*b))
        return;

    // Case 1: x is empty
    if (!body && !NW)
    {
        body = b;
        mass = b->mass;
        com_x = b->position[0];
        com_y = b->position[1];
    }
    else
    {
        // Case 2 or 3: already has body or children
        if (!NW)
        {
            // Case 3: external node then create children
            NW = new BHTree(region.NW());
            NE = new BHTree(region.NE());
            SW = new BHTree(region.SW());
            SE = new BHTree(region.SE());

            Body *existing = body;
            body = nullptr;
            insert(existing);
        }
        // still case 2: internal node then recurse
        if (NW->region.contains(*b))
            NW->insert(b);
        else if (NE->region.contains(*b))
            NE->insert(b);
        else if (SW->region.contains(*b))
            SW->insert(b);
        else if (SE->region.contains(*b))
            SE->insert(b);

        // Update center of mass and mass
        mass += b->mass;
        com_x = (com_x * (mass - b->mass) + b->mass * b->position[0]) / mass;
        com_y = (com_y * (mass - b->mass) + b->mass * b->position[1]) / mass;
    }
}

void BHTree::update_force(Body *b, double theta, double G)
{
    if (!body && !NW)
        return;
    if (body == b)
        return;

    double dx = com_x - b->position[0];
    double dy = com_y - b->position[1];
    double dist = std::sqrt(dx * dx + dy * dy + 1e-6);

    if (!NW || (region.length / dist < theta))
    {
        double F = G * b->mass * mass / (dist * dist + 1e-6);
        b->velocity[0] += F * dx / dist / b->mass;
        b->velocity[1] += F * dy / dist / b->mass;
    }
    else
    {
        NW->update_force(b, theta, G);
        NE->update_force(b, theta, G);
        SW->update_force(b, theta, G);
        SE->update_force(b, theta, G);
    }
}
