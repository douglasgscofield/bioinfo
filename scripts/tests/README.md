Tests
=====

A few simple tests for some (unfortunately not yet all) of the scripts in the
above directory.  I've haphazardly used a few hacked unit test structures and
haven't settled on any single one, but the one in intervalBed_test.sh is
probably the easiest to use.  To run tests, simply do

    make

In the top-level directory of this repository you can do the same to run these tests.

    ( cd $(git rev-parse --show-toplevel) && make )

