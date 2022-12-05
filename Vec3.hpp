// Compute the determinant of a 3x3 matrix
inline double det(double a, double b, double c, double d, double e, double f, double g, double h, double i) {
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
}
struct Vec3 {
    double x, y, z;

    // Default constructor
    Vec3() : x(0), y(0), z(0) {}

    // Constructor
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    // Overload the () operator to construct a Vec3 from an array
    Vec3(const double arr[3]) : x(arr[0]), y(arr[1]), z(arr[2]) {}

    // Overload the () operator to construct a Vec3 from a vector
    Vec3(const std::vector<double>& vec) : x(vec[0]), y(vec[1]), z(vec[2]) {}

    // Overload the () operator to construct a Vec3 from a scalar
    Vec3(double s) : x(s), y(s), z(s) {}

    // Overload the + operator to perform vector addition
    Vec3 operator+(const Vec3& other) const { return Vec3(x + other.x, y + other.y, z + other.z); }

    // Overload the - operator to perform vector subtraction
    Vec3 operator-(const Vec3& other) const { return Vec3(x - other.x, y - other.y, z - other.z); }

    // Overload the * operator to perform dot product
    double operator*(const Vec3& other) const { return x * other.x + y * other.y + z * other.z;}

    // Overload the / operator to perform cross product
    Vec3 operator/(const Vec3& other) const { return Vec3(y * other.z - z * other.y,z * other.x - x * other.z,x * other.y - y * other.x);}

    // Overload the * operator to perform scalar multiplication
    Vec3 operator*(double s) const { return Vec3(x * s, y * s, z * s); }

    // Overload the / operator to perform scalar division
    Vec3 operator/(double s) const { return Vec3(x / s, y / s, z / s); }

    // Overload the - operator to perform component-wise negation
    Vec3 operator-() const { return Vec3(-x, -y, -z); }

    // Overload the *= operator to perform vector multiplication and assignment
    Vec3& operator*=(const Vec3& other) { x *= other.x; y *= other.y; z *= other.z; return *this; }

    // Overload the /= operator to perform vector division and assignment
    Vec3& operator/=(const Vec3& other) { x /= other.x; y /= other.y; z /= other.z; return *this; }

    // Overload the += operator to perform vector addition and assignment
    Vec3& operator+=(const Vec3& other) { x += other.x; y += other.y; z += other.z; return *this; }

    // Overload the -= operator to perform vector subtraction and assignment
    Vec3& operator-=(const Vec3& other) { x -= other.x; y -= other.y; z -= other.z; return *this; }

    // Overload the -= operator for Vec3 and scalar objects
    Vec3& operator-=(double s) { x -= s; y -= s; z -= s; return *this; }

    // Overload the *= operator for Vec3 and scalar objects
    Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }

    // Overload the /= operator for Vec3 and scalar objects
    Vec3& operator/=(double s) { x /= s; y /= s; z /= s; return *this; }

    // Overload the = operator for Vec3 and scalar objects
    Vec3& operator=(double s) { x = s; y = s; z = s; return *this; }

    // Overload the [] operator to access vector elements
    double operator[](int i) const { return (*this)(i); }

    // Overload the [] operator to modify vector elements
    double& operator[](int i) { return (*this)(i); }

    // Overload the () operator to access vector elements
    double operator()(int i) const {
        if (i == 0) return x;
        if (i == 1) return y;
        if (i == 2) return z;
        throw std::out_of_range("Vec3 index out of range");
    }

    // Overload the () operator to modify vector elements
    double& operator()(int i) {
        if (i == 0) return x;
        if (i == 1) return y;
        if (i == 2) return z;
        throw std::out_of_range("Vec3 index out of range");
    }


    // Compute the matrix product of a Vec3 object and a 3x3 matrix
    friend Vec3 operator*(const Vec3& v, const double mat[3][3]) {
        return Vec3(dot(v, Vec3(mat[0][0], mat[1][0], mat[2][0])),
            dot(v, Vec3(mat[0][1], mat[1][1], mat[2][1])),
            dot(v, Vec3(mat[0][2], mat[1][2], mat[2][2])));
    }
 
    // Compute the eigenvalues and eigenvectors of a 3x3 matrix
    friend void eigen(const double mat[3][3], Vec3& eigval, Vec3 eigvec[3]) {
        // Compute the coefficients of the characteristic polynomial of mat
        double a = -mat[0][0] - mat[1][1] - mat[2][2];
        double b = mat[0][0] * mat[1][1] + mat[1][1] * mat[2][2] + mat[2][2] * mat[0][0]
            - mat[0][1] * mat[1][0] - mat[1][2] * mat[2][1] - mat[2][0] * mat[0][2];
        double c = det(mat[0][0], mat[0][1], mat[0][2],
            mat[1][0], mat[1][1], mat[1][2],
            mat[2][0], mat[2][1], mat[2][2]);
        double d = -det(mat[0][0], mat[0][1], mat[0][2],
            mat[1][0], mat[1][1], mat[1][2],
            mat[2][0], mat[2][1], mat[2][2]);

        // Compute the roots of the characteristic polynomial using the cubic formula
        double p = (3 * a * c - b * b) / (3 * a * a);
        double q = (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / (27 * a * a * a);
        double r = q / 2;
        double s = std::sqrt(r * r + (p / 3) * (p / 3) * (p / 3));
        std::complex<double> t = std::pow(std::abs(r + s), 1.0 / 3.0);
        std::complex<double> u = std::pow(std::abs(r - s), 1.0 / 3.0);
        std::complex<double> x1 = -1.0 / 3.0 * (t + u);
        // Compute the other two roots by taking the conjugate of x1
        std::complex<double> x2 = -1.0 / 3.0 * (t + u * std::complex<double>(0, 1));
        std::complex<double> x3 = -1.0 / 3.0 * (t + u * std::complex<double>(0, -1));

        // Store the eigenvalues in eigval
        eigval.x = x1.real();
        eigval.y = x2.real();
        eigval.z = x3.real();

        // Compute the eigenvectors
        for (int i = 0; i < 3; i++) {
            double mat1[3][3] = {
                {mat[0][0] - eigval.x, mat[0][1], mat[0][2]},
                {mat[1][0], mat[1][1] - eigval.x, mat[1][2]},
                {mat[2][0], mat[2][1], mat[2][2] - eigval.x}
            };
            double mat2[3][3] = {
                {mat[0][0] - eigval.y, mat[0][1], mat[0][2]},
                {mat[1][0], mat[1][1] - eigval.y, mat[1][2]},
                {mat[2][0], mat[2][1], mat[2][2] - eigval.y}
            };
            double mat3[3][3] = {
                {mat[0][0] - eigval.z, mat[0][1], mat[0][2]},
                {mat[1][0], mat[1][1] - eigval.z, mat[1][2]},
                {mat[2][0], mat[2][1], mat[2][2] - eigval.z}
            };

            // Compute the determinants of the matrices
            double det1 = det(mat1[0][0], mat1[0][1], mat1[0][2],
                mat1[1][0], mat1[1][1], mat1[1][2],
                mat1[2][0], mat1[2][1], mat1[2][2]);
            double det2 = det(mat2[0][0], mat2[0][1], mat2[0][2],
                mat2[1][0], mat2[1][1], mat2[1][2],
                mat2[2][0], mat2[2][1], mat2[2][2]);
            double det3 = det(mat3[0][0], mat3[0][1], mat3[0][2],
                mat3[1][0], mat3[1][1], mat3[1][2],
                mat3[2][0], mat3[2][1], mat3[2][2]);

            // Store the eigenvectors in eigvec
            eigvec[i].x = det1;
            eigvec[i].y = det2;
            eigvec[i].z = det3;
            eigvec[i].normalize();
        }
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

 
    // Overload the sin function for Vec3 objects
    Vec3 sin(const Vec3& v) {
        return Vec3(std::sin(v.x), std::sin(v.y), std::sin(v.z));
    }

    // Overload the cos function for Vec3 objects
    Vec3 cos(const Vec3& v) {
        return Vec3(std::cos(v.x), std::cos(v.y), std::cos(v.z));
    }

    // Overload the tan function for Vec3 objects
    Vec3 tan(const Vec3& v) {
        return Vec3(std::tan(v.x), std::tan(v.y), std::tan(v.z));
    }

    // Overload the exp function for Vec3 objects
    Vec3 exp(const Vec3& v) {
        return Vec3(std::exp(v.x), std::exp(v.y), std::exp(v.z));
    }

    // Overload the log function for Vec3 objects
    Vec3 log(const Vec3& v) {
        return Vec3(std::log(v.x), std::log(v.y), std::log(v.z));
    }

    // Overload the sqrt function for Vec3 objects
    Vec3 sqrt(const Vec3& v) {
        return Vec3(std::sqrt(v.x), std::sqrt(v.y), std::sqrt(v.z));
    }

    // Overload the pow function for Vec3 objects
    Vec3 pow(const Vec3& v, double x) {
        return Vec3(std::pow(v.x, x), std::pow(v.y, x), std::pow(v.z, x));
    }

    // Overload the abs function for Vec3 objects
    Vec3 abs(const Vec3& v) {
        return Vec3(std::abs(v.x), std::abs(v.y), std::abs(v.z));
    }

    // Overload the max function for Vec3 objects
    Vec3 max(const Vec3& v1, const Vec3& v2) {
        return Vec3(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z));
    }

    // Overload the min function for Vec3 objects
    Vec3 min(const Vec3& v1, const Vec3& v2) {
        return Vec3(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z));
    }

    // Compute the square root of a vector
    Vec3 sqrt() const {
        return Vec3(std::sqrt(x), std::sqrt(y), std::sqrt(z));
    }

    // Compute the power function of a vector and a scalar
    Vec3 pow(double s) const {
        return Vec3(std::pow(x, s), std::pow(y, s), std::pow(z, s));
    }

    // Compute the trigonometric functions of a vector
    Vec3 sin() const {
        return Vec3(std::sin(x), std::sin(y), std::sin(z));
    }
    Vec3 cos() const {
        return Vec3(std::cos(x), std::cos(y), std::cos(z));
    }
    Vec3 tan() const {
        return Vec3(std::tan(x), std::tan(y), std::tan(z));
    }

    // Compute the hyperbolic functions of a vector
    Vec3 sinh() const {
        return Vec3(std::sinh(x), std::sinh(y), std::sinh(z));
    }
    Vec3 cosh() const {
        return Vec3(std::cosh(x), std::cosh(y), std::cosh(z));
    }
    Vec3 tanh() const {
        return Vec3(std::tanh(x), std::tanh(y), std::tanh(z));
    }

    // Compute the determinant of a matrix of Vec3 objects
    friend double det(const Vec3& a, const Vec3& b, const Vec3& c) {
        return a.x * (b.y * c.z - b.z * c.y)
            - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }

    // Compute the dot product of two Vec3 objects
    friend double dot(const Vec3& a, const Vec3& b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    // Compute the cross product of two Vec3 objects
    friend Vec3 cross(const Vec3& a, const Vec3& b) {
        return Vec3(a.y * b.z - a.z * b.y,
            a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x);
    }
    // Print the vector to standard output
    void print() const {
        std::cout << "(" << x << ", " << y << ", " << z << ")" << std::endl;
    }
 
};
