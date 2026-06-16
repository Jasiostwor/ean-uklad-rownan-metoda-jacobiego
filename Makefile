CXX = clang++
CXXFLAGS = -std=c++17 -Wall `wx-config --cxxflags` -I/opt/homebrew/include
LDFLAGS = `wx-config --libs` -L/opt/homebrew/lib -lmpfr -lgmp

TARGET = JacobiApp
SRCS = main.cpp Jacobi.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
