# refactor

The goal of refactoring is:
* Add helper functions or methods as needed
* Simplify logic
* Assess the structure of the code for possible architectural issues
* Assess the design pattern and propose new changes if needed
* Prompt me on how to proceed if multiple design patterns may fit

Do not perform tiny changes. The point is to address big picture issues.

Helper functions should not be created if they will be 2 lines or less.

Helper functions should ideally take primitives for parameters and return types.
Primitives include built-in Python data types.
