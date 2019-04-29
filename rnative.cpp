#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector> 

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

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

DLLIntializer init_hook{};


extern "C"
SEXP basicFunction(SEXP x, SEXP y)
{
	std::cerr << " [TRACE] Calling R from C++" << std::endl;
	// Allocate a vector of size 1 of type real (double)
	// Remember: in R, a scalar is a vector of size one  
	SEXP out = PROTECT(allocVector(REALSXP, 1));
	// asReal(x) => Interprets the first element of the
	// vector x as a double 
	REAL(out)[0] = 10.0 * asReal(x) + 4.0 * asReal(y);
	//     UNPROTECT(< Number of allocations> )
	// OR: UNPROTECT(< Number of protects> )
	UNPROTECT(1);	
	return out;
}

extern "C"
SEXP GetVersion()
{
	SEXP version;
	// Allocate 1 string 
	PROTECT(version = allocVector(STRSXP, 1));
	SET_STRING_ELT(version, 0, ::mkChar("version 0.1"));
	UNPROTECT(1);
	return version;
}

extern "C"
auto metrics(SEXP xs) -> SEXP 
{
	size_t size = ::length(xs);
	double* vec =  REAL(xs);
	double sumsq  = 0.0;
	double sum    = 0.0;
	double x;
	for(auto i = 0; i < size; i++){
		x = vec[i];
		sumsq = sumsq + x * x;
		sum   = sum + x;
	}
	std::cout << std::setprecision(4) << std::fixed;
	std::cout << std::endl;
	std::cout << " =>    Mean = " << sum / size << std::endl;
	std::cout << " =>    Norm = " << std::sqrt(sumsq) << std::endl;
	std::cout << " =>    Size = " << size << std::endl;
	// If the function returns nothing, it should return R_NilValue 
	return R_NilValue;
}

/** Computes the equation: 
  * For each i = 0, 1, 2, ... N - 1
  *   zs_vector[i] = a * xs_vector[i] + b * ys_vector[i] + c 
  *
  * Returns: Vector zs_vector[i]
  */
extern "C"
SEXP linearcomb(SEXP xs_vector, SEXP ys_vector, SEXP a, SEXP b, SEXP c)
{
	size_t size = ::length(xs_vector);
	// Allocate output vector of same size 
	SEXP output = PROTECT(allocVector(REALSXP, size));
	double *pout, *pxs, *pys;
	pout = REAL(output);
	pxs  = REAL(xs_vector);
	pys  = REAL(ys_vector);
	for(size_t i = 0; i < size; i++)
		pout[i] = asReal(a) * pxs[i] + asReal(b) * pys[i] + asReal(c);
	UNPROTECT(1);
	return output;
}

#if 0 
extern "C"
SEXP computeStatistics(SEXP xs_vector)
{
	size_t size = ::length(xs_vector);
	
	SEXP mean = PROTECT(allocVector(REALSXP, 1));
	SEXP norm = PROTECT(allocVector(REALSXP, 1));	
	
	double sum, sumsq, x;
	sum = 0.0, sumsq = 0.0, x = 0.0;
	for(size_t i = 0; i < size; i++){
		x = REAL(xs_vector)[i];
		sum += x;
		sumsq += x * x;
	}
	REAL(mean)[0] = sum / size;
	REAL(norm)[0] = std::sqrt(sumsq);

	SEXP list_names = PROTECT(allocVector(STRSXP, 2));
	// Set list names 
	SET_STRING_ELT(list_names, 0, mkChar("mean"));
	SET_STRING_ELT(list_names, 1, mkChar("norm"));

	// Set list values (vectors)
	SEXP list = PROTECT(allocVector(VECSXP, 2));	
	SET_VECTOR_ELT(list, 0, mean);
	SET_VECTOR_ELT(list, 1, norm);	
	setAttrib(list, R_NameSymbol, list_names);

	// UNPROTECT( <<Number of PROTECTs or of allocVector>> )
	UNPROTECT(4);	
	return list;
}
#endif 
