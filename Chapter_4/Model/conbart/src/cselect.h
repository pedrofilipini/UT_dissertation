#ifndef c_select_h
#define c_select_h

#include "arma_config.h"
#include <RcppArmadillo.h>

#include <vector>
#include "funs.h"
#include "info.h"
#include <tgmath.h>

//Here is a nice solution provided by
//Luc Touraille on Stack Overflow
//(https://stackoverflow.com/questions/698520/search-for-nearest-value-in-an-array-of-doubles-in-c/701141#701141)

template <typename BidirectionalIterator, typename T>
BidirectionalIterator getClosest(BidirectionalIterator first,
                                 BidirectionalIterator last,
                                 const T & value)
{
  BidirectionalIterator before = std::lower_bound(first, last, value);

  if (before == first) return first;
  if (before == last)  return --last; // iterator must be bidirectional

  BidirectionalIterator after = before;
  --before;

  return (*after - value) < (value - *before) ? after : before;
}

template <typename BidirectionalIterator, typename T>
std::size_t getClosestIndex(BidirectionalIterator first,
                            BidirectionalIterator last,
                            const T & value)
{
  return std::distance(first, getClosest(first, last, value));
}
//End of the Stack Overflow solution


class c_select {
public:
  virtual double val(double x) = 0;
};

class c_select_interp: public c_select {
  public:
  std::vector<double> data;
  std::vector<double> c_data;
  size_t data_length;

  c_select_interp(std::vector<double> data_, std::vector<double> c_data_, size_t data_length_) {
    data = data_;
    c_data = c_data_;
    data_length = data_length_;
  }
  c_select_interp() {
    data_length = 0;
  }
  double val(double x) {
    int idx;
    double result;
    //Search closest value in data for a given x
    //Save the index of that value
    //Return y of that index

    idx = getClosestIndex(&data[0], &data[data_length], x);

    //Rcpp::Rcout << "idx " << idx << endl;

    //Add the interpolation
    if(x<=data[0]){ //Case extreme left
      result = c_data[0];
    }else if(x>=data[(data_length-1)]){ //case extreme right
      //Make Linear Regression
      double b =(c_data[(data_length-1)]-c_data[(data_length-2)])/(data[(data_length-1)]-data[(data_length-2)]);
      double a = c_data[(data_length-2)]-(b*data[(data_length-2)]);
      result = a+b*x;
      //Rcpp::Rcout << "a: " << a << "; b: " << b << endl;
    }else if(x<=data[idx]){ //middle case where x is closest to upper idx
      //Make Linear Regression
      double b =(c_data[idx]-c_data[(idx-1)])/(data[idx]-data[(idx-1)]);
      double a = c_data[(idx-1)]-(b*data[(idx-1)]);
      result = a+b*x;
    }else if(x>data[idx]){ //middle case where x is closest to lower idx
      //Make Linear Regression
      double b =(c_data[idx+1]-c_data[(idx)])/(data[idx+1]-data[(idx)]);
      double a = c_data[(idx)]-(b*data[(idx)]);
      result = a+b*x;
    }

    return(result);
  }
};

#endif
