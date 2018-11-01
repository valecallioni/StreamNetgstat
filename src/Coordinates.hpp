#ifndef PACSPROJECT_COORDINATES_HPP
#define PACSPROJECT_COORDINATES_HPP

class Coordinates{

private:
    double x1;
    double x2;

public:
    Coordinates() = default;
    Coordinates(double x, double y):x1(x), x2(y){};
    const double getX1() const {return x1;};
    const double getX2() const {return x2;};
};

#endif //PACSPROJECT_COORDINATES_HPP
