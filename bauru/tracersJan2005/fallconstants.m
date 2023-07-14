function [ai,bi,gi,as,bs,gs,ag,bg,gg]=fallconstants(constants)

%alternative ice fall speed constants
switch constants(1)
case 1
    ai=71.3;
    bi=0.6635;
    gi=0.5;
case 2
    ai=700;
    bi=1;
    gi=0.35;
case 3
    ai=3.25;
    bi=0.33;
    gi=0.3;
end

switch constants(2)
case 1
    as=4.84;
    bs=0.25;
    gs=0.5;
case 2
    as=8.97;
    bs=0.42;
    gs=0.5;
case 3
    as=1.14;
    bs=0.11;
    gs=0.4;
case 4
    as=27.7;
    bs=0.5;
    gs=0.5;
case 5
    as=17.0;
    bs=0.5;
    gs=0.5;
case 6
    as=6.96;
    bs=0.33;
    gs=0.3;
end


switch constants(3)
case 1
    ag=253;
    bg=0.734;
    gg=0.422;    
case 2
    ag=94.5;
    bg=0.5;
    gg=0.5;
case 3
    ag=19.3;
    bg=0.37;
    gg=0.5;
case 4
    ag=124;
    bg=0.64;
    gg=0.5;
case 5
    ag=199;
    bg=0.8;
    gg=0.286;
    
end



    
    