#pragma once
#include "body.h"
#include <gtkmm.h>
#include <vector>
#include "barnes_hutt.h"

class Visualizer : public Gtk::Window {
public:
    Visualizer(const double* all_positions, int n_bodies, int steps, const std::vector<BHTree*>& trees);

protected:
    Gtk::DrawingArea drawing_area;
    const double* all_positions;
    int n_bodies;
    int total_frames;
    int current_frame = 0;
    const BHTree* tree = nullptr;
    std::vector<BHTree*> trees;

    bool on_draw(const Cairo::RefPtr<Cairo::Context>& cr) override;
    bool on_timeout();
};

