CXX = g++
CXXFLAGS = `pkg-config --cflags --libs gtkmm-4.0` -std=c++17 -pthread -Iinclude 
SRC = ssrc/main.cpp src/body.cpp src/visualizer.cpp src/simulation.cpp
OBJ = obj/main.o obj/body.o
TARGET = n_body_simulation

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

obj/%.o: src/%.cpp | obj
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj:
	mkdir -p obj

clean:
	rm -rf obj $(TARGET)