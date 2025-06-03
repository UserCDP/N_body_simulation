#pragma once
#include <gtkmm.h>
// #include "body.h"
#include <vector>
#include <mutex>

class Visualizer : public Gtk::Window {
public:
    Visualizer(std::vector<Body>* bodies, std::mutex* mutex);

protected:
    bool on_draw(const Cairo::RefPtr<Cairo::Context>& cr) override;
    bool on_timeout();

    Gtk::DrawingArea drawing_area;
    std::vector<Body>* bodies;
    std::mutex* mutex;
};