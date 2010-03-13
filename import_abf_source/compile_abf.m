% script to compile and link abf interface MEX code
% JAB 6/15/07w

% cd 'C:\Program Files\MATLAB\R2006a\work\Nanopore\Current\import_abf\import_abf_source\build'
!del import_abf.mexw32
%mex -D_abf_DEBUG_ -D_abf_USING_MEX_ -I. -L. import_abf.c abf_interface.c abffio_omf.lib
mex -D_abf_USING_MEX_ -I. -L. import_abf.c abf_interface.c abffio_omf.lib
% mex -D_abf_USING_MEX_ -I. -L. import_abf.cpp abf_interface.cpp abffio_omf.lib

%!cp 'C:\Program Files\MATLAB\R2006a\work\Nanopore\Jan10_2008\import_abf\import_abf_source\build\import_abf.mexw32' 'C:\Program Files\MATLAB\R2006a\work\Nanopore\Jan10_2008\import_abf\import_abf.mexw32'
