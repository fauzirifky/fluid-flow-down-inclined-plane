// make_tuple example
#include <iostream>
#include <tuple>
#include <functional>
#include <string>
#include <array>
using namespace std;

std::tuple<double, double> fun(double x, double ){
  double a = 0;
  for (const auto & xj : x ){
    a = a + xj;
  }
  return make_tuple(a, a);
}



int main(){
  double balance[] = {1.2, 1.3, 1.1, 1.6, 20.0};
  // Declaring tuple
  //auto first = std::make_tuple (10,'a');             // tuple < int, char >

  //const int a = 0; int b[3];                         // decayed types:
  //auto second = std::make_tuple (a,b);               // tuple < int, int* >

  //auto third = std::make_tuple (std::ref(a),"abc");  // tuple < const int&, const char* >

  //std::cout << "third contains: " << std::get<0>(third);
  //std::cout << " and " << std::get<1>(third);
  //std::cout << std::endl;
  auto v = fun(balance);
  double num1 = std::get<0>(v);
  double num2 = std::get<1>(v);
  std::cout<< "size of the array is " <<num2<<", with average "<<num1<<"\n";
}
