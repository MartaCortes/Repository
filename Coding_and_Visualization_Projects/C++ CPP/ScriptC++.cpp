#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

double my_knn_C2_inverse_euclidean(NumericMatrix X,NumericVector X0,NumericVector y){
  
  int nrows = X.nrow();
  
  
  double closest_distance1 = 99999999;
  double closest_output1 = -1;
  int closest_neighbor1 = -1;
  
  for(int i=0; i<nrows; ++i){
    double distance = 0;
    distance = sum(pow((X(i,_) - X0), 2.0));
    distance = sqrt(distance);
    
    if(distance < closest_distance1){
      closest_distance1 = distance;
      closest_output1 = y[i];
      closest_neighbor1 = i;
    }
  }
  
  double closest_distance2 = 99999999;
  double closest_output2 = -1;
  int closest_neighbor2 = -1;
  
  for(int i=0; i<nrows; ++i){
    double distance = 0;
    distance = sum(pow((X(i,_) - X0), 2.0));
    distance = sqrt(distance);
    
    if((distance < closest_distance2) && (distance!=closest_distance1)){
      closest_distance2 = distance;
      closest_output2 = y[i];
      closest_neighbor2 = i;
    }
  }
  
  NumericVector closest_output = NumericVector::create(closest_output1,closest_output2);
  double result_closest_output = (closest_output1/closest_distance1+closest_output2/closest_distance2)/((1/closest_distance1)+(1/closest_distance2));
  return result_closest_output;

  
}  