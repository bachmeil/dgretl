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
	dmd test4.d ../gretl/*.d ../gretlexp/dynmodel.d -L-lgretl-1.0 -version=standalone
	./test4
	rm test4

5logit:
	dmd logit.d ../gretl/*.d -L/usr/local/lib/libgretl-1.0.so -version=standalone
	./logit
	rm logit

6probit:
	dmd probit.d ../gretl/*.d -L/usr/local/lib/libgretl-1.0.so -version=standalone
	./probit
	rm probit

7adf:
	dmd adf.d ../gretl/*.d ../gretlexp/dynmodel.d -L-lgretl-1.0 -version=standalone
	./adf
	rm adf

dmvnorm:
	dmd dmvnorm.d ../gretl/*.d -L-lgretl-1.0 -version=standalone
	./dmvnorm
	rm dmvnorm
