function d=fdlatlon

d=fzero(@zero,-45.554);


function x=zero(lat);

lat1=-19.015;
lon1=-52.498;

x=distlatlon(lat1,lon1,lat,lon1)-sqrt(1000^2/2);

