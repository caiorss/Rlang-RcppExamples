#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <string>
#include <map>

#include <Rcpp.h>


struct DLLIntializer{
	DLLIntializer()
	{
		std::cerr << "[TRACE] - DLL loaded into process OK." << std::endl;
	}
	~DLLIntializer()
	{
		std::cerr << "[TRACE] - DLL unloaded from process OK." << std::endl;
	}
};

// Global object allocated in process' data segment for
// DLL logging during loading and unloading. 
DLLIntializer init_hook{};

RcppExport
auto libraryVersion() -> SEXP
{
    return Rcpp::wrap("version 1.0 - alpha");
}

/** Compute equation: 
 * for_each(i) => ys[i] = a * xs[i] + b 
 * return ys 
 */
RcppExport
SEXP linearcomb(SEXP a, SEXP b, SEXP xs)
{
	Rcpp::NumericVector vx(xs);    
    double aa = Rcpp::as<double>(a);
    double bb = Rcpp::as<double>(b);
	Rcpp::NumericVector out(vx.size());	
	size_t n = vx.size();
	// std::cerr << " [TRACE] Size of xs = " << n << std::endl;
	for(auto i = 0; i < vx.size(); i++)
		out[i] = aa * vx[i] + bb;
	return out;
}

RcppExport
auto computeStatistics(SEXP xs) -> SEXP
{
    Rcpp::NumericVector vx{xs};
    double sum = 0.0;
    double sumsq = 0.0;

    std::for_each(vx.begin(), vx.end(), [&](double x)
    {
        sum   += x;
        sumsq += x * x;
    });

    auto result = std::map<std::string, double>
    {
        {"sum",   sum},
        {"sumsq", sumsq},
        {"mean",  sum / vx.size()},
        {"norm",  std::sqrt(sumsq)}
    };
    return Rcpp::wrap(result);
}

RcppExport
auto computeStatistics2(SEXP xs) -> SEXP
{
    Rcpp::NumericVector vx{xs};
    double sum = 0.0;
    double sumsq = 0.0;

    std::for_each(vx.begin(), vx.end(), [&](double x)
    {
        sum   += x;
        sumsq += x * x;
    });
    return Rcpp::List::create(
                    Rcpp::Named("sum", sum)
                   ,Rcpp::Named("sumsq", sumsq)
                   ,Rcpp::Named("mean", sum / vx.size())
                   ,Rcpp::Named("norm", std::sqrt(sumsq))
                );
}


RcppExport
SEXP tabulate(SEXP mfunc, SEXP mvec)
{
    Rcpp::Function     func = mfunc;
    Rcpp::NumericVector vec = mvec;
    std::cout << std::setprecision(4) << std::fixed;
    for(auto const& x: vec)
        std::cout << std::setw(10) << x
                  << std::setw(10) << Rcpp::as<double>(func(x))
                  << std::endl;
    // Return nothing on R side
    return R_NilValue;
}


/* Compute:
 *
 *  FOR i = 0 TO N DO:
 *    out[i] = sqrt(va[i] * 4.0 + vb[i] * 5.0
 *  RETURN ouy
 *
 * Note: The vectorized operation are not efficient
 * Rcpp::sqrt(4.0 * va + 5.0 * vb) is equivaent to:
 * Rcpp::sqrt(operatior*(4.0, va) + operator*(5.0, vb))
 * Rcpp::sqrt(operator+(operatior*(4.0, va), operator*(5.0, vb)))
 *
 * Every operator function call is 1 loop, so, the total number
 * of passes over the input data is equal to:
 *  =>  1 (SQRT Call) + 3 (Number of operators) = 4
 *
 * This function requires 4 for-loops over the data,
 * thus a single for-loop is more efficient.
 *
 */
RcppExport
SEXP vectorOperations(SEXP sa, SEXP sb)
{
    using Vector = Rcpp::NumericVector;
    Vector va = sa;
    Vector vb = sb;
    Vector result = Rcpp::sqrt(4.0 * va + 5.0 * vb);
    return result;
}

RcppExport
SEXP ShowNormalRandoms(SEXP s_size)
{
    size_t n = Rcpp::as<size_t>(s_size);
    std::cout << " [INFO] size = " << n << std::endl;

    // Get R function rnorm for generating N normally distributed
    // random numbers
    Rcpp::Function rnorm("rnorm");

    // Call R built-in function rnorm
    // that generates n normally distributed random numbers with
    // mean 20.0 and standard deviation std = 5.0
    //---------------------------------------------
    // NOTE: Every R fuction returns an S-expression struct
    SEXP sexp = rnorm(Rcpp::Named("n", n),
                      Rcpp::Named("sd", 5.0),
                      Rcpp::Named("mean", 20.0)
                      );
    // Convert s-expression to numeric vector
    Rcpp::NumericVector result{sexp};

    std::cout << std::setprecision(5) << std::fixed;
    size_t k = 0;
    std::for_each(result.begin(), result.end(),
                  [&](auto x) -> void
                   {
                    std::cout << std::setw(10) << k++
                              << std::setw(10) << x
                              << std::endl;
                   });
    // std::for_each(result.begin())
    return R_NilValue;
}




