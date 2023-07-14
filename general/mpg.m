%works out the mpg given the miles and the number of litres
function mpg=mpg(miles,litres)

L2g=1/0.586/8; %conversion factor for litres to gallons

mpg=miles/(L2g*litres);