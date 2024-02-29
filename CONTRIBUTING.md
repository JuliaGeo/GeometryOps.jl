# Contributing to GeometryOps.jl

Thank you for your interest in contributing to GeometryOps.jl! Here are some guidelines to help you get started.

## Code Organization

- Try to keep each functionality in one file or folder.  For example, polygon cutting and intersection are in different files.  However, common utilities were moved out to another file, which is also OK.  Semantically, each file should have one job.
- We use a [**literate programming**](https://en.wikipedia.org/wiki/Literate_programming) approach.  This means that the code and documentation are intertwined, making the code easier to understand. Please follow this style when contributing.
    - Please add an example in the comments at the top of each file. This should illustrate how to use the functionality provided by the file.  We've found that these examples are incredibly helpful for understanding the code.
    - Visual examples using plots (we use Makie in this library) are always helpful, and should go at the top of the file.  Each file should have a 
- If you are confused, take a look at e.g. `src/methods/distance.jl` or `src/methods/area.jl`, and follow their structure as templates.

## Tests

Tests are great!  We love tests.  They are organized in the `test` directory in the same structure as the `src` directory.  Try and add as many tests as possible, and try for integration tests - especially those with real-world uses!


