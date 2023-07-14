ichoose_styles=1;

clear line_pattern

for idat2=1:length(labs)
    switch labs(idat2).l
        case {'LWP day AMSRE time3 0pt03'}
            line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
            labs(idat2).l = '3%';
        case {'LWP day AMSRE time3 0pt05'}
            line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 1 0]; marker_style(istyle).m='d'; istyle=istyle+1;
            labs(idat2).l = '5%';            
        case {'LWP day AMSRE time3 0pt07'}
            line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='d'; istyle=istyle+1;
            labs(idat2).l = '7%';
        case {'LWP day AMSRE time3 0pt09'}
            line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
            labs(idat2).l = '9%';            
        case {'LWP day AMSRE time3 0pt2'}
            line_pattern(istyle).p= '--';  line_colour(istyle).c=[1 0.65 0]; marker_style(istyle).m='d'; istyle=istyle+1;
            labs(idat2).l = '20%';            
                             
        
            
    end
end