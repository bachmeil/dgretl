1:
	dmd test1.d ../gretl/*.d -L/usr/local/lib/libgretl-1.0.so -version=standalone
	./test1
	rm test1

2:
	dmd test2.d ../gretl/*.d -L/usr/local/lib/libgretl-1.0.so -version=standalone
	./test2
	rm test2

3ts:
	dmd test3.d ../gretl/*.d ../gretlexp/dynmodel.d -L/usr/local/lib/libgretl-1.0.so -version=standalone
	./test3
	rm test3

4ts:
	dmd test4.d ../gretl/*.d ../gretlexp/dynmodel.d -L/usr/local/lib/libgretl-1.0.so -version=standalone
	./test4
	rm test4
