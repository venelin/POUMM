Fixed the old date and missing reference in Description field. 

This is a minor bug-fix, which should ensure compatibility with data.table v1.12.2 (to be released soon). 

R CMD CHECK passes with 0 Errors, 0 Warnings and 1 Note:

checking installed package size ... NOTE
  installed size is  6.1Mb
  sub-directories of 1Mb or more:
    libs   4.6Mb

This is due to the relatively big size of the compiled C++ code for the package. I think this size would vary, depending on the C++ compiler used. There is nothing I can do to reduce it.



