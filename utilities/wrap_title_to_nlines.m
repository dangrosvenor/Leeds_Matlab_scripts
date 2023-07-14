function [tit_wrapped] = wrap_title_to_nlines(tit_str,ncols,nlines)

%ncols=55; nlines=3; %nlines is the no. of lines we are aiming for (will add blank ones to make up to this later)
tit_wrapped = textwrap({tit_str},ncols);

%Add extra lines if they don't exist
for iw=length(tit_wrapped)+1:nlines
    tit_wrapped{iw}='';
end
title(tit_wrapped);

