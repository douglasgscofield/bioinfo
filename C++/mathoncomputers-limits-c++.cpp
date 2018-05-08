#include <iostream>
#include <limits>

#define intsize(_t_) std::cout << "size of " #_t_ " is " << sizeof(_t_) << " bytes"  \
    << " min " << +std::numeric_limits<_t_>::min()                                   \
    << " max " << +std::numeric_limits<_t_>::max()                                   \
    << std::endl;

int main() {
    intsize(signed char);  // unary '+' coaxes std::cout to print chars as numeric
    intsize(short int);
    intsize(int);
    intsize(long int);
    intsize(long long int);
    intsize(unsigned char);
    intsize(unsigned short int);
    intsize(unsigned int);
    intsize(unsigned long int);
    intsize(unsigned long long int);
}
