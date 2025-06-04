#include "visualizer.h"
#include <cmath>
#include <iostream>
#include <mutex>


Visualizer::Visualizer(const double* all_positions, int n_bodies, int total_frames)
    : all_positions(all_positions), n_bodies(n_bodies), total_frames(total_frames)
{
    set_title("N-Body Simulation");
    set_default_size(800, 800);
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

bool Visualizer::on_draw(const Cairo::RefPtr<Cairo::Context>& cr) {
    const double scale = 1.0 / 1e5;  
    const double offset = 400;
    const double radius = 3.0;

    cr->set_source_rgb(0, 0, 0);
    cr->paint();

    cr->set_source_rgb(1, 1, 1);

    for (int i = 0; i < n_bodies; ++i) {
        double x = all_positions[(current_frame * n_bodies + i) * 2] * scale + offset;
        double y = all_positions[(current_frame * n_bodies + i) * 2 + 1] * scale + offset;

        if (x >= 0 && y >= 0 && x <= 800 && y <= 800) {
            cr->arc(x, y, radius, 0, 2 * M_PI);
            cr->fill();
        }
    }
    
    ++current_frame;
    return true;
}


