# Overview

This module provides an interface to the [Gretl](http://gretl.sourceforge.net/) API.

What functionality does it provide? Among other things, it offers:

- Matrix operations (multiplication, transpose, determinant, etc.)
- Descriptive statistics
- Random number generation and distribution functions
- Regression analysis, including generalized linear models

It does a lot more than that. These are the parts I use.

The original version of this library was written in 2013, and since it
has done what I needed it to do, other libraries may now be available
that do things better. (I don't know; I haven't had a reason to check.)

# Very Brief History

When I started using D in 2013, I couldn't find anything that fit my
needs, so I started creating bindings to libgretl in my spare time. I
chose Gretl because it provided the most convenient interface to LAPACK.
Moreover, it was created as an econometrics library, so the more specialized
tools fit with what I need to do. I wrote the bindings, added some
operator overloading to get convenient syntax for the basic stuff, and
I was good to go.

# Approach

I've found that it helps to add a section like this to a public project.
I'm a believer that the only sensible way to write this type of code is
to use the garbage collector. In the cases where that leads to noticeable
slowdowns, you can always handle the memory management manually using the
functionality provided by Gretl.[^1] In short, my philosophy is that the
garbage collector is your friend. Some minor modifications of the code
is necessary if you disagree.

[^1]: Note that Gretl is a C library, so all memory management is done
manually. The functions to manually allocate and free data structures is
available if you want to go that route. Alternatively, it's easy to do
reference counting with D - and I've done that elsewhere - but it's not
clear there will be any efficiency benefits to doing so. RC adds a lot of
overhead. My experience suggests that the real gains come from avoiding
the allocation of memory rather than how you manage the freeing of the
memory. Pauses in your program are a complete non-issue for these applications.
We're not writing video games.

There are two ways to make use of a fast language like D. First, you can
use it to write highly optimized code that pushes hardware performance to
its limits. Second, you can use it to write code in a fun, convenient
manner that runs slower than the optimized stuff, but more than fast
enough. 

My work falls into the latter camp because my objective function
puts weight on developer time and the time it takes to iterate. Simple,
specialized functions are preferred to generic code most of the time.
The underlying libraries need to be easy to read and modify. Testing is
important, but it's even more important to write code that's easy to test
and less likely to have bugs.

# Basic Usage

This example (coming soon!) works on Ubuntu Linux 18.04. There's no reason it won't work
with proper modifications elsewhere.

Although I will use Dub, since that's what the kids use to build D projects
these days, I don't currently have plans to put it on code.dlang.org.
It's open source, so someone willing to volunteer the time should do so.
You'll have to do it the old fashioned way, by cloning this repo and
telling Dub where to find it. See below (coming soon!) for Makefile usage.

# Individual Modules

- [base.d](doc-base.html)
