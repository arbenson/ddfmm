#include "file_io.h"

#include <stdio>
#include <vector>

int main() {
#if 0
    std::vector<Point3> points;
    std::vector<Point3> coords;
	
    ReadWrl("F16.wrl", points, coords);
    std::cout << "size of points (should be 1165): " << points.size() << std::endl;
    std::cout << "size of coords (should be 2162): " << coords.size() << std::endl;

    for (int i = 0; i < points.size(); ++i) {
	std::cout << points[i] << std::endl;
    }

    points.clear();
    coords.clear();

    ReadWrl("SubmarineJ.wrl", points, coords);
    std::cout << "size of points (should be 459): " << points.size() << std::endl;
    std::cout << "size of coords (should be 640): " << coords.size() << std::endl;

    points.clear();
    coords.clear();

    ReadWrl("sphere.wrl", points, coords);
    std::cout << "size of points (should be 1681): " << points.size() << std::endl;
    std::cout << "size of coords (should be 3200): " << coords.size() << std::endl;
#endif

    std::string geom = "sphere.wrl";
    double K = 1024;
    double NPW = 10;
    int NCPU = 8;
    int NC = 8;
    NewData(geom, K, NPW, NCPU, NC);
}
