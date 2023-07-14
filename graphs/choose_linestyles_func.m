function [colour,patt,mark] = choose_linestyles_func(j,cdan,pdan,markers)

size(cdan);
scdan=ans(2);

smark=length(markers);
spdan=length(pdan);

            if rem(j,scdan)==0
                colour.c=cdan(scdan).c;
            else
                colour.c=cdan(rem(j,scdan)).c;
            end
            if rem(j,spdan)==0
                patt.p=pdan(spdan).p;
            else
                patt.p=pdan(rem(j,spdan)).p;
            end
            if rem(j,smark)==0
                mark.m=markers(smark).m;
            else
                mark.m=markers(rem(j,smark)).m;
            end