# Random notes on gretl.matrix

- The module publicly imports gretl.matfunctions and gretl.vector since
    they're so closely related.

## Versions

- It is designed to work with R code. You can call Gretl functions and
    call the resulting functions from R. Use `version=r` when compiling.
- You should use `version=standalone` if you're creating a D executable
    that will call Gretl functions.
- You should use `version=inline` if you're writing functions to be
    called from R. `version=standalone` and `version=inline` conflict.

## DoubleMatrix

- This is the main matrix type. D allocates the memory, and it's managed
    by the D garbage collector. It's fully compatible with Gretl's
    matrix type, so a `DoubleMatrix` can be passed to any Gretl function
    that takes a matrix.
- If you want to avoid the garbage collector, use the `GretlMatrix`
    struct instead. You'll need to call `.free()` when you're done with
    it. Allocation of the struct (a trivial part of the memory allocation)
    is done by the garbage collector. The data itself is allocated using
    malloc. `GretlMatrix` can be passed to any Gretl function that takes
    a matrix.
- If you want a pointer to the underlying data array, that's done by 
    calling `.ptr`. If you want a pointer to a `GretlMatrix` struct that
    can be passed to a Gretl function, use `.matptr`.
    
    Example:
    
    ```
    // x is a DoubleMatrix holding data
    foo(x); // foo can take either DoubleMatrix or GretlMatrix, do to alias this
    foo(x.ptr); // pass double * to function foo
    foo(x.matptr); // pass GretlMatrix * to function foo
    foo(&x); // pass DoubleMatrix * to function foo, different from x.matptr and rarely used
    
    GretlMatrix gm;
    gm.free(); // Frees the data in the array, but the struct itself still exists
    ```
- See examples/doublemat.d for example usage of the `DoubleMatrix` struct.
    That file serves as the documentation. If there is something not covered
    in that file or if the file doesn't run properly, please create a
    bug report so that it can be fixed.
