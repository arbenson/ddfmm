#include "file_io.h"
#include "vec3t.h"

#include <stdio>
#include <vector>

#define MAX_LINE_LENGTH 256;

int ReadWrl(std::string fname, std::vector<Point3>& points,
            std::vector<Point3>& coords) {
    char filename[MAX_FILE_NAME_LENGTH];
    sprintf(filename, "data/%s", geom_file.c_str());
    FILE *geom_file;
    geom_file = fopen (filename, "r");
    if (geom_file == NULL) {
	std::cerr << "Failed to open geometry file: " << filename << std::endl;
	return -1;
    }

    points.clear();
    coords.clear();
    char pts[MAX_LINE_LENGTH]
    fscanf(geom_file, "%s", &pts);
    if (strncmp(pts, "points", 6)) {
	std::cerr << "Did not detect 'points' at the beginning of geometry file" << std::endl;
	return -1;
    }

    // Format of points is:
    //         x1   x2   x3
    double x1, x2, x3;
    while (fscanf (geom_file, "%f %f %f", &x1, &x2, &x3) == 3) {
	points.push_back(Point3(x1, x2, x3));
	std::cout << x1 << x2 << x3 << std::endl;
    }
    // Should fail on reading "coords"

    // Format of coords is:
    //          xcoord ycoord zcoord -1
    int c1, c2, c3, c4;
    while (fscanf (geom_file, "%d %d %d %d", &c1, &c2, &c3, &c4) == 4) {
	if (c4 != -1) {
	    std::cerr << "Expected -1 in fourth coordinate" << std::endl;
	    return -1;	    
	}
	coords.push_back(Point3(c1 + 1, c2 + 1, c3 + 1));
	std::cout << c1 << c2 << c3 << std::endl;
    }

    // Normalize points to [-1, 1]
    // Find min and max points along each row
    double first = points[0][0];
    Point3 min(first, first, first);
    Point3 max(first, first, first);
    for (int i = 0; i < points.size(); ++i) {
	Point3 p = points[i];
	for (int j = 0; j < 3; ++j) {
	    min[j] = std::min(min[j], p[j]);
	    max[j] = std::max(max[j], p[j]);
	}
    }
    Point3 mid = (max + min) / 2;
    Point3 diff = (max - min) / 2;
    double avg_diff = std::max(std::max(diff[0], diff[1]), diff[2]);
    for (int i = 0; i < points.size(); ++i) {
	points[i] -= mid;
	points[i] /= max_avg_diff;
    }

    return 0;
}
