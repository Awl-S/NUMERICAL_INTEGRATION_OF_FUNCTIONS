#include "src/NumericalIntegrator.h"

int main() {
    std::function<double(double)> f = [](double x) {
        return x * x; // Define the function to integrate here
    };

    RectangleMethod rectangleMethod(f);
    TrapezoidMethod trapezoidMethod(f);
    SimpsonMethod simpsonMethod(f);
    GaussMethod gaussMethod(f);

    double a = 0.0; // Lower bound of integration
    double b = 1.0; // Upper bound of integration
    int n = 1000; // Number of intervals

    std::cout << "Rectangle Method: " << rectangleMethod.Integrate(a, b, n) << std::endl;
    std::cout << "Trapezoid Method: " << trapezoidMethod.Integrate(a, b, n) << std::endl;
    std::cout << "Simpson Method: " << simpsonMethod.Integrate(a, b, n) << std::endl;
    std::cout << "Gauss Method: " << gaussMethod.Integrate(a, b, n) << std::endl;

    return 0;
}