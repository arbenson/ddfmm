#include "file_io.h"

#include <stdio>
#include <vector>

int main() {
    std::vector<Point3> points;
    std::vector<Point3> coords;
	
#if 0
    ReadWrl("F16.wrl", points, coords);
    std::cout << "size of points (should be 1165): " << points.size() << std::endl;
    std::cout << "size of coords (should be 2162): " << coords.size() << std::endl;

    for (int i = 0; i < points.size(); ++i) {
	std::cout << points[i] << std::endl;
    }

    points.clear();
    coords.clear();
#endif


    ReadWrl("SubmarineJ.wrl", points, coords);
    std::cout << "size of points (should be 459): " << points.size() << std::endl;
    std::cout << "size of coords (should be 640): " << coords.size() << std::endl;

#if 0
    points.clear();
    coords.clear();

    ReadWrl("sphere.wrl", points, coords);
    std::cout << "size of points (should be 1681): " << points.size() << std::endl;
    std::cout << "size of coords (should be 3200): " << coords.size() << std::endl;
#endif
}
