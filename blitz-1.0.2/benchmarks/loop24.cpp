
// loop24 generated by makeloops.py Thu Jun 30 16:44:56 2011

#include <blitz/vector2.h>
#include <blitz/array.h>
#include <random/uniform.h>
#include <blitz/benchext.h>

#ifdef BZ_HAVE_VALARRAY
 #define BENCHMARK_VALARRAY
#endif

#ifdef BENCHMARK_VALARRAY
#include <valarray>
#endif

namespace blitz {
extern void sink();
}

using namespace blitz;
using namespace std;

#if defined(BZ_FORTRAN_SYMBOLS_WITH_TRAILING_UNDERSCORES)
 #define loop24_f77 loop24_f77_
 #define loop24_f77overhead loop24_f77overhead_
 #define loop24_f90 loop24_f90_
 #define loop24_f90overhead loop24_f90overhead_
#elif defined(BZ_FORTRAN_SYMBOLS_WITH_DOUBLE_TRAILING_UNDERSCORES)
 #define loop24_f77 loop24_f77__
 #define loop24_f77overhead loop24_f77overhead__
 #define loop24_f90 loop24_f90__
 #define loop24_f90overhead loop24_f90overhead__
#elif defined(BZ_FORTRAN_SYMBOLS_CAPS)
 #define loop24_f77 LOOP24_F77
 #define loop24_f77overhead LOOP24_F77OVERHEAD
 #define loop24_f90 LOOP24_F90
 #define loop24_f90overhead LOOP24_F90OVERHEAD
#endif

extern "C" {
  void loop24_f77(const int& N, double* x, double* a, double* b, double* c, double* d, double* y);
  void loop24_f77overhead(const int& N, double* x, double* a, double* b, double* c, double* d, double* y);
  void loop24_f90(const int& N, double* x, double* a, double* b, double* c, double* d, double* y);
  void loop24_f90overhead(const int& N, double* x, double* a, double* b, double* c, double* d, double* y);

}

void VectorVersion(BenchmarkExt<int>& bench);
void ArrayVersion(BenchmarkExt<int>& bench);
void ArrayVersion_unaligned(BenchmarkExt<int>& bench);
void ArrayVersion_misaligned(BenchmarkExt<int>& bench);
void ArrayVersion_index(BenchmarkExt<int>& bench);
void doTinyVectorVersion(BenchmarkExt<int>& bench);
void F77Version(BenchmarkExt<int>& bench);
#ifdef FORTRAN_90
void F90Version(BenchmarkExt<int>& bench);
#endif
#ifdef BENCHMARK_VALARRAY
void ValarrayVersion(BenchmarkExt<int>& bench);
#endif

const int numSizes = 80;
const bool runvector=false; // no point as long as Vector is Array<1>

int main()
{
    int numBenchmarks = 5;
    if (runvector) numBenchmarks++;
#ifdef BENCHMARK_VALARRAY
    numBenchmarks++;
#endif
#ifdef FORTRAN_90
    numBenchmarks++;
#endif

    BenchmarkExt<int> bench("loop24: $x = $a*$c - $b*$c; $y = $a*$d + $b+$c", numBenchmarks);

    bench.setNumParameters(numSizes);

    Array<int,1> parameters(numSizes);
    Array<long,1> iters(numSizes);
    Array<double,1> flops(numSizes);

    parameters=pow(pow(2.,0.25),tensor::i)+tensor::i;
    flops = 6 * parameters;
    iters = 100000000L / flops;
    iters = where(iters<2, 2, iters);
    cout << iters << endl;
    
    bench.setParameterVector(parameters);
    bench.setIterations(iters);
    bench.setOpsPerIteration(flops);
    bench.setDependentVariable("flops");
    bench.beginBenchmarking();



    ArrayVersion(bench);
    ArrayVersion_unaligned(bench);
    ArrayVersion_misaligned(bench);
    ArrayVersion_index(bench);
    //doTinyVectorVersion(bench);
    F77Version(bench);
#ifdef FORTRAN_90
    F90Version(bench);
#endif
#ifdef BENCHMARK_VALARRAY
    ValarrayVersion(bench);
#endif

    if(runvector)
      VectorVersion(bench);

    bench.endBenchmarking();

    bench.saveMatlabGraph("loop24.m");
    return 0;
}

template<class T>
void initializeRandomDouble(T* data, int numElements, int stride = 1)
{
    ranlib::Uniform<T> rnd;

    for (int i=0; i < numElements; ++i)
        data[size_t(i*stride)] = rnd.random();
}

template<class T>
void initializeRandomDouble(valarray<T>& data, int numElements, int stride = 1)
{
    ranlib::Uniform<T> rnd;

    for (int i=0; i < numElements; ++i)
        data[size_t(i*stride)] = rnd.random();
}

void VectorVersion(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Vector<T>");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();
        long iters = bench.getIterations();

        cout << bench.currentImplementation() << ": N = " << N << endl;

        Vector<double> x(N);
        initializeRandomDouble(x.data(), N);
        Vector<double> a(N);
        initializeRandomDouble(a.data(), N);
        Vector<double> b(N);
        initializeRandomDouble(b.data(), N);
        Vector<double> c(N);
        initializeRandomDouble(c.data(), N);
        Vector<double> d(N);
        initializeRandomDouble(d.data(), N);
        Vector<double> y(N);
        initializeRandomDouble(y.data(), N);


        bench.start();
        for (long i=0; i < iters; ++i)
        {
            x = a*c - b*c; y = a*d + b+c;
            sink();
        }
        bench.stop();

        bench.startOverhead();
        for (long i=0; i < iters; ++i) {
            sink();
	}

        bench.stopOverhead();
    }

    bench.endImplementation();
}


  void ArrayVersion(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Array<T,1>");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();
        long iters = bench.getIterations();

        cout << bench.currentImplementation() << ": N = " << N << endl;

        Array<double,1> x(N);
        initializeRandomDouble(x.dataFirst(), N);
        Array<double,1> a(N);
        initializeRandomDouble(a.dataFirst(), N);
        Array<double,1> b(N);
        initializeRandomDouble(b.dataFirst(), N);
        Array<double,1> c(N);
        initializeRandomDouble(c.dataFirst(), N);
        Array<double,1> d(N);
        initializeRandomDouble(d.dataFirst(), N);
        Array<double,1> y(N);
        initializeRandomDouble(y.dataFirst(), N);


        bench.start();
        for (long i=0; i < iters; ++i)
        {
            x = a*c - b*c; y = a*d + b+c;
            sink();
        }
        bench.stop();

        bench.startOverhead();
        for (long i=0; i < iters; ++i) {
            sink();
	}

        bench.stopOverhead();
    }

    bench.endImplementation();
}


  void ArrayVersion_index(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Array<T,1> (indexexpr.)");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();
        long iters = bench.getIterations();

        cout << bench.currentImplementation() << ": N = " << N << endl;

        Array<double,1> x(N);
        initializeRandomDouble(x.dataFirst(), N);
        Array<double,1> a(N);
        initializeRandomDouble(a.dataFirst(), N);
        Array<double,1> b(N);
        initializeRandomDouble(b.dataFirst(), N);
        Array<double,1> c(N);
        initializeRandomDouble(c.dataFirst(), N);
        Array<double,1> d(N);
        initializeRandomDouble(d.dataFirst(), N);
        Array<double,1> y(N);
        initializeRandomDouble(y.dataFirst(), N);


        bench.start();
        for (long i=0; i < iters; ++i)
        {
            x = a(tensor::i)*c(tensor::i) - b(tensor::i)*c(tensor::i); y = a(tensor::i)*d(tensor::i) + b(tensor::i)+c(tensor::i);;
            sink();
        }
        bench.stop();

        bench.startOverhead();
        for (long i=0; i < iters; ++i) {
            sink();
	}

        bench.stopOverhead();
    }

    bench.endImplementation();
}

  void ArrayVersion_unaligned(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Array<T,1> (unal.)");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();
        long iters = bench.getIterations();

        cout << bench.currentImplementation() << ": N = " << N << endl;


    Array<double,1> xfill(N+1);
    Array<double,1> x(xfill(Range(1,N)));
    initializeRandomDouble(x.dataFirst(), N);

    Array<double,1> afill(N+1);
    Array<double,1> a(afill(Range(1,N)));
    initializeRandomDouble(a.dataFirst(), N);

    Array<double,1> bfill(N+1);
    Array<double,1> b(bfill(Range(1,N)));
    initializeRandomDouble(b.dataFirst(), N);

    Array<double,1> cfill(N+1);
    Array<double,1> c(cfill(Range(1,N)));
    initializeRandomDouble(c.dataFirst(), N);

    Array<double,1> dfill(N+1);
    Array<double,1> d(dfill(Range(1,N)));
    initializeRandomDouble(d.dataFirst(), N);

    Array<double,1> yfill(N+1);
    Array<double,1> y(yfill(Range(1,N)));
    initializeRandomDouble(y.dataFirst(), N);


        bench.start();
        for (long i=0; i < iters; ++i)
        {
            x = a*c - b*c; y = a*d + b+c;
            sink();
        }
        bench.stop();

        bench.startOverhead();
        for (long i=0; i < iters; ++i) {
            sink();
	}

        bench.stopOverhead();
    }

    bench.endImplementation();
}

  void ArrayVersion_misaligned(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Array<T,1> (misal.)");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();
        long iters = bench.getIterations();

        cout << bench.currentImplementation() << ": N = " << N << endl;


    Array<double,1> xfill(N+6);
    Array<double,1> x(xfill(Range(0,N+0-1)));
    initializeRandomDouble(x.dataFirst(), N);

    Array<double,1> afill(N+6);
    Array<double,1> a(afill(Range(1,N+1-1)));
    initializeRandomDouble(a.dataFirst(), N);

    Array<double,1> bfill(N+6);
    Array<double,1> b(bfill(Range(2,N+2-1)));
    initializeRandomDouble(b.dataFirst(), N);

    Array<double,1> cfill(N+6);
    Array<double,1> c(cfill(Range(3,N+3-1)));
    initializeRandomDouble(c.dataFirst(), N);

    Array<double,1> dfill(N+6);
    Array<double,1> d(dfill(Range(4,N+4-1)));
    initializeRandomDouble(d.dataFirst(), N);

    Array<double,1> yfill(N+6);
    Array<double,1> y(yfill(Range(5,N+5-1)));
    initializeRandomDouble(y.dataFirst(), N);


        bench.start();
        for (long i=0; i < iters; ++i)
        {
            x = a*c - b*c; y = a*d + b+c;
            sink();
        }
        bench.stop();

        bench.startOverhead();
        for (long i=0; i < iters; ++i) {
            sink();
	}

        bench.stopOverhead();
    }

    bench.endImplementation();
}

#ifdef BENCHMARK_VALARRAY
void ValarrayVersion(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("valarray<T>");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();
        cout << bench.currentImplementation() << ": N = " << N << endl;

        long iters = bench.getIterations();

        valarray<double> x(N);
        initializeRandomDouble(x, N);
        valarray<double> a(N);
        initializeRandomDouble(a, N);
        valarray<double> b(N);
        initializeRandomDouble(b, N);
        valarray<double> c(N);
        initializeRandomDouble(c, N);
        valarray<double> d(N);
        initializeRandomDouble(d, N);
        valarray<double> y(N);
        initializeRandomDouble(y, N);


        bench.start();
        for (long i=0; i < iters; ++i)
        {
            x = a*c - b*c; y = a*d + b+c;
            sink();
        }
        bench.stop();

        bench.startOverhead();
        for (long i=0; i < iters; ++i) {
	  sink();
	}
        bench.stopOverhead();
    }

    bench.endImplementation();
}
#endif

void F77Version(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Fortran 77");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();
        cout << bench.currentImplementation() << ": N = " << N << endl;

        int iters = bench.getIterations();

        double* x = new double[N];
        initializeRandomDouble(x, N);
        double* a = new double[N];
        initializeRandomDouble(a, N);
        double* b = new double[N];
        initializeRandomDouble(b, N);
        double* c = new double[N];
        initializeRandomDouble(c, N);
        double* d = new double[N];
        initializeRandomDouble(d, N);
        double* y = new double[N];
        initializeRandomDouble(y, N);
        

        bench.start();
        for (int iter=0; iter < iters; ++iter)
            loop24_f77(N, x, a, b, c, d, y);
        bench.stop();

        bench.startOverhead();
        for (int iter=0; iter < iters; ++iter)
            loop24_f77overhead(N, x, a, b, c, d, y);

        bench.stopOverhead();

        delete [] x;
        delete [] a;
        delete [] b;
        delete [] c;
        delete [] d;
        delete [] y;

    }

    bench.endImplementation();
}

#ifdef FORTRAN_90
void F90Version(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Fortran 90");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();
        cout << bench.currentImplementation() << ": N = " << N << endl;

        int iters = bench.getIterations();

        double* x = new double[N];
        initializeRandomDouble(x, N);
        double* a = new double[N];
        initializeRandomDouble(a, N);
        double* b = new double[N];
        initializeRandomDouble(b, N);
        double* c = new double[N];
        initializeRandomDouble(c, N);
        double* d = new double[N];
        initializeRandomDouble(d, N);
        double* y = new double[N];
        initializeRandomDouble(y, N);


        bench.start();
        for (int iter=0; iter < iters; ++iter)
            loop24_f90(N, x, a, b, c, d, y);
        bench.stop();

        bench.startOverhead();
        for (int iter=0; iter < iters; ++iter)
            loop24_f90overhead(N, x, a, b, c, d, y);

        bench.stopOverhead();
        delete [] x;
        delete [] a;
        delete [] b;
        delete [] c;
        delete [] d;
        delete [] y;

    }

    bench.endImplementation();
}
#endif

