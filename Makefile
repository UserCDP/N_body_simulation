CXX = g++
CXXFLAGS = -std=c++11 -pthread -Iinclude
SRC = src/main.cpp src/body.cpp
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