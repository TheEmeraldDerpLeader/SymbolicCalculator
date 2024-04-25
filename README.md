# SymbolicCalculator
Small C++ project to represent multivariable polynomials and perform computations on them.

Represent polynomials (i.e. symbolic expressions) with SymExp objects, which are a vector of Product objects. Each Product object contains ids and powers to efficiently represent factors of the product. 

#Current Features:
- multiplying, adding, and subtraction expressions with each other. 

- Substituting variables with other SymExps using the SymExpTable class, which can be used to evaluate the symbolic expressions like a function. Also supports printing out expressions by converting symbol ids to text, converting back into SymExp is planned, but low priority.

- Calculating the derivative or gradient of expressions

- Using generalized Newton's method to find zeros of a expression or vector of expressions (to represent a system of polynomials) with quadratic convergence near roots

- Displaying the fractals that result from Newton's method using SFML and textures

##Future Features

My current focus is on flattening SymExp into a class that stores all data in a single region of memory. Currently SymExp uses vectors of vectors, throwing in another vector when you need multiple SymExps to represent a system. All of this means a substantial amount of computation time for evaluation is spent on memory calls rather than computation, which moving everything into a single dynamic array should hopefully reduce, allowing for faster generation of fractals. After this, I plan to port the code over to OpenCL, allowing for SymExp evaluation on GPUs. To make this simple, the FlatSymExp class will be designed in a way that would easily port over to C.