CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra
TARGET   = quadratura
SRC      = main.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)

.PHONY: all run clean
