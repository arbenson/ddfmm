#include "file_io.h"

#include <stdio>
#include <vector>

int main() {
    std::vector<Point3> points;
    std::vector<Point3> coords;

    ReadWrl("F16.wrl", points, coords);
    for (int i = 0; i < points.size(); ++i) {
	std::cout << points[i] << std::endl;
    }
    for (int i = 0; i < coords.size(); ++i) {
	std::cout << coords[i] << std::endl;
    }
}
