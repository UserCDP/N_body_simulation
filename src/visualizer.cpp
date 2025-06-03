#include "visualizer.h"
#include <cmath>
#include <iostream>
#include <mutex>

/**
 * Constructor for the Visualizer class for the graphical representation of the simulation.
 * @param b     The list of bodies in the simulation
 * @param m     A mutex
 */
Visualizer::Visualizer(std::vector<Body>* b, std::mutex* m)
    : bodies(b), mutex(m) {
    set_title("N-Body Simulation");
    set_default_size(800, 800);
    add(drawing_area);
    drawing_area.signal_draw().connect(sigc::mem_fun(*this, &Visualizer::on_draw));
    drawing_area.show();
    Glib::signal_timeout().connect(sigc::mem_fun(*this, &Visualizer::on_timeout), 33);
}

bool Visualizer::on_timeout() {
    drawing_area.queue_draw();
    return true;
}

bool Visualizer::on_draw(const Cairo::RefPtr<Cairo::Context>& cr) {
    const double radius = 3.0;
    const double scale = 6.0;
    cr->set_source_rgb(0, 0, 0);
    cr->paint();

    std::lock_guard<std::mutex> lock(*mutex);
    cr->set_source_rgb(1, 1, 1);
    for (const auto& body : *bodies) {
        double x = body.position[0] * scale;
        double y = body.position[1] * scale;
        cr->arc(x, y, radius, 0, 2 * M_PI);
        cr->fill();
    }
    static int frame = 0;
    if (frame < 300) {
        auto surface = cr->get_target();
        surface->flush();

        Cairo::RefPtr<Cairo::ImageSurface> img = Cairo::RefPtr<Cairo::ImageSurface>::cast_dynamic(surface);
        std::string tempname = "frames/tmp_" + std::to_string(frame) + ".png";
        std::string finalname = "frames/frame_" + std::to_string(frame) + ".png";

        try {
            img->write_to_png(tempname);
            std::rename(tempname.c_str(), finalname.c_str());
            frame++;
        } catch (...) {
            std::cerr << "Skipped frame " << frame << " due to PNG write error\n";
        }
    }
    return true;
}
