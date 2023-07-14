year=0;
loan=14188;

while loan>0
    year=year+1/12;
    loan=loan + 103 - 4.8/12*loan;
end

year