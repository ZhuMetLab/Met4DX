
Dynamically linked libraries compiled for Windows (32-bit and 64-bit), and Linux (64-bit)
are contained in this archive.

To use the Windows libraries successfully, you need to have the "Visual C++
Redistributable for Visual Studio 2017" installed on your system.

* schema-version_3_0.sql - documentation of the database schema defined in analysis.tdf files

* include/c/timsdata.h - definition and documentation of the C API.

* include/c/tsfdata.h - definition and documentation of the C API for spectrum data.

* examples/py contains an example python wrapper and a small example program.
  To run these examples, you will need a python intallation with numpy and mathplot.

* examples/timsdataSampleCpp contains a C++ example (for Visual Studio 2017).
  Please also note the readme file in thirdparty if you try to compile the
  C++ example program. If you want to run the C++ example program on windows 
  copy timsdata.dll (win32) to the folder of timsdataSampleCpp.exe.

The file "redist.txt" contains a list of files that may be distributed together with
software using tdf-sdk. It also lists files that must be distributed in such cases.

Further documentation on specific fields or properties in the database are provided upon request.

