% m-file to import some les parameters
FileName=input('Enter the file name','s');
i=input('Enter the array number');
cd ../LESIMPORT

% run the les import progam 2-d only
Importdiag3dnew
cd ../LESTOOLS
% run the les processing program uses i
AFTERIMPORTDIAG2;

% array of usefull temperatures
V=[0 -5 -8 -15 -20 -28 -35 -40 -50 -60 -70];


