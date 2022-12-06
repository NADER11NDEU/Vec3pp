# Vec3pp
A Vec3 object represents a 3-dimensional vector with x, y, and z components. The Vec3 class provides a number of methods for performing operations on 3-dimensional vectors, such as vector addition, subtraction, dot product, cross product, and scalar multiplication and division. Vec3 objects can be constructed from arrays, vectors, scalars, or individual x, y, and z components. The class also overloads a number of operators, such as +, -, *, and /, to allow for intuitive and concise vector operations.

# Features
* Overloaded operators for vector operations such as addition, subtraction, dot product, cross product, and scalar multiplication and division.
* Functions for computing the magnitude and normalizing the vector.
* Function for computing the eigenvalues and eigenvectors of a 3x3 matrix.

# Usage
To use the Vec3 class, include the Vec3.hpp header file in your project and create an instance of the class using the Vec3 constructor. The class provides several overloads for the +, -, *, and / operators, as well as assignment operators and other utility functions.

# Example
To use the Vec3 class, simply include the Vec3.hpp header file and create a Vec3 object with the appropriate coordinates:

```c++
#include "Vec3.hpp"

Vec3 v(1, 2, 3);
``` 

You can then use the overloaded operators to perform vector operations:


```c++
Vec3 a(1, 2, 3);
Vec3 b(4, 5, 6);

// Vector addition
Vec3 c = a + b;

// Vector subtraction
Vec3 d = a - b;

// Dot product
double e = a * b;

// Cross product
Vec3 f = a / b;

// Scalar multiplication
Vec3 g = a * 2.0;

// Scalar division
Vec3 h = a / 2.0;

```

The Vec3 class also provides functions for computing the magnitude and normalizing the vector:
```c++
Vec3 v(1, 2, 3);

// Compute the magnitude of the vector
double magnitude = v.mag();

// Normalize the vector
Vec3 normalized = v.norm();
```

