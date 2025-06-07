CXX = g++
CXXFLAGS = -std=c++17 `pkg-config --cflags --libs gtkmm-3.0` -pthread -Iinclude 
SRC = src/main.cpp src/body.cpp src/visualizer.cpp src/simulation.cpp src/barnes_hutt.cpp
OBJ = obj/main.o obj/body.o obj/visualizer.o obj/simulation.o obj/barnes_hutt.o
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