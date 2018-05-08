Documentation practices group
=============================

I'm attempting to put this file in the `$HOME/Box Sync` folder on my desktop, to see if it gets synced to my laptop and especially back to the Documentation practices group folder on Box.

It does... but the Box Sync client on my laptop had to be restarted to get the new folder synced there.  My limited experience so far suggests it is not as robust as the Dropbox client is in this regard.

## Syncing files to your desktop/laptop

Files will not automatically be synced to your desktop.  The `Documentation practices group` folder was *not* synced to my desktop/laptop until I did the following:

1. In the box website (<https://uppsala.app.box.com>), at All Files, hover to the right size of the folder name until the three-dots "More Options" box appears.  Click there.
2. A dropdown will appear.  Click More Actions.
3. A popup will appear.  Select Sync.

This should be enough to get the file in the list of "Synced" files on the left of the All Files window and thus one of the files Box thinks it should keep synced.

This is also the **only** way to get files Synced to your desktop from the Box site.  Much better to do this to an entire folder, so its contents will be tracked.

## Code formatting

I agree with Linnea that the lack of attractive code formatting is a bummer:

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

Compare this to a copy of this README on Github: <https://github.com/douglasgscofield/bioinfo/tree/master/C++>

and the C++ file itself: <https://github.com/douglasgscofield/bioinfo/blob/master/C++/mathoncomputers-limits-c++.cpp>

I've also placed a copy of the C++ file within the Box folder to view it on the Box website.  It has a little extra formatting, but it is not very well done.
