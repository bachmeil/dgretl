library(embedr)
compileFile("sampling.d", "foo", "dmdgretl")
.Call("genobs", 10L)
.Call("genobs", 1000L)
