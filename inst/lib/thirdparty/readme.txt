For the C++ sample, you need:
- boost (http://www.boost.org/)
- CppSqLite: https://github.com/neosmart/CppSQLite
- SqLite itself (https://sqlite.org/download.html, use sqlite-amalgamation-<version>.zip)

Download the archives from these the sites mentioned above, copy them into this directory and unpack them.
Rename the directory sqlite-amalgamation-<version> to sqlite.
Rename the directory boost_<version> to boost.

You should get the following directory structure 
 - CppSQLite-master
   - (three files in the directory)
 - sqlite
   - (four files in the directory)
 - boost
   - boost
   - libs
   - (much more)
 
You do not have to compile any of these third party projects for compiling the sample program.
However when writing a "real" program using timsdata, you might choose to compile some or all of
these projects separately.