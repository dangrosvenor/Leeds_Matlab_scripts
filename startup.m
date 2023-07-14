%script which is executed on startup of Matlab.
%set(0,'Defaultpaperpositionmode','auto');
%above doesn't work as only applies to figures

set(0,'defaultaxesyminortick','on');
set(0,'defaultaxesxminortick','on'); %makes Matlab plot small ticks in between the large ones

%for netCDF SNCtools
%javaaddpath ( 'C:\Program Files\MATLAB\R2007b\toolbox\toolsUI-2.2.22.jar' );
%setpref ( 'SNCTOOLS', 'USE_JAVA', true );