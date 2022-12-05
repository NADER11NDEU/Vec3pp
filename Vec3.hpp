class Vec3 {
 
public:
    double x;
    double y;
    double z;
    // Constructor
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    // Overload the + operator to perform vector addition
    Vec3 operator+(const Vec3& other) const {
        return Vec3(x + other.x, y + other.y, z + other.z);
    }

    // Overload the - operator to perform vector subtraction
    Vec3 operator-(const Vec3& other) const {
        return Vec3(x - other.x, y - other.y, z - other.z);
    }

    // Overload the * operator to perform dot product
    double operator*(const Vec3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    // Overload the / operator to perform cross product
    Vec3 operator/(const Vec3& other) const {
        return Vec3(y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x);
    }

    // Overload the * operator to perform scalar multiplication
    Vec3 operator*(double s) const {
        return Vec3(x * s, y * s, z * s);
    }

    // Overload the / operator to perform scalar division
    Vec3 operator/(double s) const {
        return Vec3(x / s, y / s, z / s);
    }

    // Overload the - operator to perform component-wise negation
    Vec3 operator-() const {
        return Vec3(-x, -y, -z);
    }

    // Overload the *= operator to perform vector multiplication and assignment
    Vec3& operator*=(const Vec3& other) {
        x *= other.x;
        y *= other.y;
        z *= other.z;
        return *this;
    }

    // Overload the /= operator to perform vector division and assignment
    Vec3& operator/=(const Vec3& other) {
        x /= other.x;
        y /= other.y;
        z /= other.z;
        return *this;
    }

    // Overload the += operator to perform vector addition and assignment
    Vec3& operator+=(const Vec3& other) {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }

    // Overload the -= operator to perform vector subtraction and assignment
    Vec3& operator-=(const Vec3& other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }
    // Return the magnitude of the vector
    double mag() const {
        return std::sqrt(mag2());
    }

    // Return the square of the magnitude of the vector
    double mag2() const {
        return x * x + y * y + z * z;
    }

    // Normalize the vector
    Vec3 normalize() const {
        double m = mag();
        return Vec3(x / m, y / m, z / m);
    }

    // Return the angle between two vectors
    double angle(const Vec3& other) const {
        return std::acos(*this * other / mag() / other.mag());
    }

    // Return the projection of this vector onto another vector
    Vec3 project(const Vec3& other) const {
        return other * (*this * other) / other.mag2();
    }

    // Return the component of this vector in the direction of another vector
    Vec3 component(const Vec3& other) const {
        return other.normalize() * (*this * other.normalize());
    }

    // Rotate the vector around the x axis by the given angle
    Vec3 rotateX(double angle) const {
        double s = std::sin(angle);
        double c = std::cos(angle);
        return Vec3(x, c * y - s * z, s * y + c * z);
    }

    // Rotate the vector around the y axis by the given angle
    Vec3 rotateY(double angle) const {
        double s = std::sin(angle);
        double c = std::cos(angle);
        return Vec3(c * x + s * z, y, -s * x + c * z);
    }

    // Rotate the vector around the z axis by the given angle
    Vec3 rotateZ(double angle) const {
        double s = std::sin(angle);
        double c = std::cos(angle);
        return Vec3(c * x - s * y, s * x + c * y, z);
    }

    // Rotate the vector around an arbitrary axis by the given angle
    Vec3 rotate(const Vec3& axis, double angle) const {
        double s = std::sin(angle);
        double c = std::cos(angle);
        return (*this * c) + (axis * (*this * axis)) * (1 - c) + (*this / axis) * s;
    }

    // Return the minimum component of the vector
    double min() const {
        return std::min(std::min(x, y), z);
    }

    // Return the maximum component of the vector
    double max() const {
        return std::max(std::max(x, y), z);
    }

    // Clamp the vector to the specified range
    Vec3 clamp(double min, double max) const {
        return Vec3(std::clamp(x, min, max),
            std::clamp(y, min, max),
            std::clamp(z, min, max));
    }

    // Compute the distance between two vectors
    double distance(const Vec3& other) const {
        return std::sqrt(distance2(other));
    }

    // Compute the square of the distance between two vectors
    double distance2(const Vec3& other) const {
        return (x - other.x) * (x - other.x) + (y - other.y) * (y - other.y) + (z - other.z) * (z - other.z);
    }

    // Check if two vectors are equal within a specified tolerance
    bool equal(const Vec3& other, double epsilon) const {
        return std::abs(x - other.x) < epsilon &&
            std::abs(y - other.y) < epsilon &&
            std::abs(z - other.z) < epsilon;
    }

    // Check if two vectors are not equal within a specified tolerance
    bool notEqual(const Vec3& other, double epsilon) const {
        return !equal(other, epsilon);
    }

    // Compute the linear interpolation between two vectors
    Vec3 lerp(const Vec3& other, double t) const {
        return (*this * (1 - t)) + (other * t);
    }

 
    // Print the vector to standard output
    void print() const {
        std::cout << "(" << x << ", " << y << ", " << z << ")" << std::endl;
    }
 
};
