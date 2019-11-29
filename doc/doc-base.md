# `gretl.base`

This module will almost always have to be imported. It includes the
infrastructure that supports everything else.

Of greatest interest are the definitions of some structs:

- `Dataset` is used to hold information about a dataset you've created or
    loaded from a file. The specific usage is determined by what you're
    doing.
- `Model` is used to hold the output of different types of model estimation.

Some functions:

- `loadDataset` is used to load a dataset from a file. *Warning:* The underlying
    Gretl function allocates a dataset and returns a pointer to it. If
    you want to load data from a file and call Gretl functions for further
    analysis, this is what you need. You'll have to call the `destroy_dataset`
    function on it when you're done. If you don't want to mess with this,
    call `loadMatrix` instead.
- `loadMatrix` copies the data from a `Dataset *` into a `DoubleMatrix`.
    It would be
    inconvenient to have to manage the pointer returned by Gretl (and
    thus you'd have to free the memory manually when you're done), so
    it copies the data into a `DoubleMatrix` struct for you, after which
    it destroys the `Dataset` allocated by Gretl.
- `seq` is a convenience function like R's seq function. The
    start and end arguments are inclusive.
    
    Example:
    
    ```
    seq(0, 6) // equivalent to [0, 1, 2, 3, 4, 5, 6]
    ```

There's also a bunch of C function prototypes. None of the other modules
would work without them.
