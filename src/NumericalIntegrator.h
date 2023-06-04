#ifndef NUMERICAL_INTEGRATION_OF_FUNCTIONS_NUMERICALINTEGRATOR_H
#define NUMERICAL_INTEGRATION_OF_FUNCTIONS_NUMERICALINTEGRATOR_H

#include <cmath>
#include <functional>
#include <vector>
#include <iostream>
#include <gsl/gsl_integration.h>

class NumericalIntegrator {
public:
    NumericalIntegrator(const std::function<double(double)>& f) : f_(f) {}
    virtual ~NumericalIntegrator() {}

    virtual double Integrate(double a, double b, int n) const = 0;

protected:
    std::function<double(double)> f_;
};

class RectangleMethod : public NumericalIntegrator {
public:
    using NumericalIntegrator::NumericalIntegrator;
    double Integrate(double a, double b, int n) const override {
        double h = (b - a) / n;
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += f_(a + h * (i + 0.5)); // midpoint
        }
        return h * sum;
    }
};

class TrapezoidMethod : public NumericalIntegrator {
public:
    using NumericalIntegrator::NumericalIntegrator;
    double Integrate(double a, double b, int n) const override {
        double h = (b - a) / n;
        double sum = 0.5 * (f_(a) + f_(b));
        for (int i = 1; i < n; i++) {
            sum += f_(a + i * h);
        }
        return h * sum;
    }
};

class SimpsonMethod : public NumericalIntegrator {
public:
    using NumericalIntegrator::NumericalIntegrator;
    double Integrate(double a, double b, int n) const override {
        double h = (b - a) / n;
        double sum = f_(a) + f_(b);
        for (int i = 1; i < n; i++) {
            sum += (i % 2 == 0 ? 2 : 4) * f_(a + i * h);
        }
        return h / 3 * sum;
    }
};

class GaussMethod {
public:
    GaussMethod(const std::function<double(double)>& f) : f_(f) {}

    double Integrate(double a, double b, int n) const {
        gsl_integration_workspace * w = gsl_integration_workspace_alloc(n);

        double result, error;
        gsl_function F;
        F.function = &GaussMethod::Wrapper;
        F.params = (void*)&f_;

        gsl_integration_qags(&F, a, b, 0, 1e-7, n, w, &result, &error);

        gsl_integration_workspace_free(w);

        return result;
    }

private:
    static double Wrapper(double x, void* params) {
        std::function<double(double)>* f = static_cast<std::function<double(double)>*>(params);
        return (*f)(x);
    }

    std::function<double(double)> f_;
};



#endif //NUMERICAL_INTEGRATION_OF_FUNCTIONS_NUMERICALINTEGRATOR_H
