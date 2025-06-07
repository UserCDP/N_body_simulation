#include "visualizer.h"
#include <cmath>
#include <iostream>
#include <mutex>
const int SIZE = 800;


Visualizer::Visualizer(const double* all_positions, int n_bodies, int total_frames, const BHTree* final_tree)
    : all_positions(all_positions), n_bodies(n_bodies), total_frames(total_frames), tree(final_tree)
{
    set_title("N-Body Simulation");
    set_default_size(SIZE, SIZE);
    add(drawing_area);
    drawing_area.signal_draw().connect(sigc::mem_fun(*this, &Visualizer::on_draw));
    drawing_area.show();
    Glib::signal_timeout().connect(sigc::mem_fun(*this, &Visualizer::on_timeout), 33);
}

bool Visualizer::on_timeout() {
    if (current_frame >= total_frames) return false;
    drawing_area.queue_draw();
    return true;
}


void draw_tree(const Cairo::RefPtr<Cairo::Context>& cr, const BHTree* node, double scale, double offset_x, double offset_y) {
    if (!node) return;
    double x = node->region.cx;
    double y = node->region.cy;
    double len = node->region.length;
    // Convert to screen coordinates
    double left   = (x - len/2) * scale + offset_x;
    double top    = (y - len/2) * scale + offset_y;

    cr->set_source_rgba(0.3, 0.6, 1.0, 0.2);
    cr->rectangle(left, top, len * scale, len * scale);
    cr->stroke();

    draw_tree(cr, node->NW, scale, offset_x, offset_y);
    draw_tree(cr, node->NE, scale, offset_x, offset_y);
    draw_tree(cr, node->SW, scale, offset_x, offset_y);
    draw_tree(cr, node->SE, scale, offset_x, offset_y);
}

bool Visualizer::on_draw(const Cairo::RefPtr<Cairo::Context>& cr) {
    const double scale = 1.0 / 1e13; 
    Gtk::Allocation alloc = get_allocation();
    const int width = alloc.get_width();
    const int height = alloc.get_height();
    const double offset_x = width / 2.0;
    const double offset_y = height / 2.0;
    const double radius = 3.0;
    
    cr->set_source_rgb(0, 0, 0);
    cr->paint();
    // Draw tree boxes
    if (tree) {
        draw_tree(cr, tree, scale, offset_x, offset_y);
    }
    // Draw bodies
    cr->set_source_rgb(1, 1, 1);

    for (int i = 0; i < n_bodies; ++i) {
        double x = all_positions[(current_frame * n_bodies + i) * 2] * scale + offset_x;
        double y = all_positions[(current_frame * n_bodies + i) * 2 + 1] * scale + offset_y;

        if (x >= 0 && y >= 0 && x <= width && y <= height) {
            cr->arc(x, y, radius, 0, 2 * M_PI);
            cr->fill();
        }
    }
    
    ++current_frame;
    return true;
}


