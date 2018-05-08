Documentation practices group
=============================

I'm attempting to put this file in the `$HOME/Box Sync` folder on my desktop, to see if it gets synced to my laptop and especially back to the Documentation practices group folder on Box.

I agree with Linnea that the lack of code formatting is a bummer:

```cpp
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
```

Compare this to a copy on Github: <https://github.com/douglasgscofield/bioinfo/blob/master/C%2B%2B/mathoncomputers-limits-c%2B%2B.cpp>

