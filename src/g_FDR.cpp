#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector g_FDR(NumericMatrix fdr){
	int nrow = fdr.nrow(), ncol = fdr.ncol();
	NumericVector temp_result(nrow * ncol);
	NumericMatrix result(nrow, ncol);	
	
	/*Convert the input matrix into a vector*/
	NumericVector temp(nrow * ncol);
	int k = 0;
	for(int r = 0; r < nrow; ++r){
		for(int c = 0; c < ncol; ++c){
			temp[k++] = fdr(r, c);
		}
	}
	int n = temp.size();
	
	/*Calculate the global FDR*/
	for(int i = 0; i < n; ++i){
		double total = 0;
		int num = 0;
		double u = temp[i];
		for(int j = 0; j < n; ++j){			
			if(temp[j] <= u){
				total += temp[j];
				num += 1;
			}			
		}
		double FDR = total / num;
		temp_result[i] = FDR;
	}
	
	/*Return the result matrix*/
	int l = 0;
	for(int r = 0; r < nrow; ++r){
		for(int c = 0; c < ncol; ++c){
			result(r, c) = temp_result[l++];
		}
	}
	return result;
}