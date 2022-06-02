#include <fftw3.h>

#include <Eigen/Core>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

using namespace std;
using namespace Eigen;

typedef Eigen::Matrix<double, 6, 1> Vec6;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Vector2d Vec2;

const double pi = M_PI;