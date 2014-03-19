#include "file_io.h"

#include <stdio>
#include <vector>

int main() {
    std::string geomfile = "F16.wrl";
    double K = 512;
    double NPW = 10;
    int NCPU = 4;
    int NC = 8;
    IntNumTns geom;

    NewData(geomfile, K, NPW, NCPU, NC, geom);

    std::cout << "Geometry-------" << std::endl
	      << geom << std::endl;
}
