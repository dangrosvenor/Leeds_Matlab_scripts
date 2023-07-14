ichoose_styles=1;

for idat2=1:length(labs)
    switch labs(idat2).l
        case {'TLWP day AMSRE','LWP day AMSRE'}
            line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
        case {'TLWP night AMSRE','LWP night AMSRE'}
            line_pattern(istyle).p= '--';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
        case {'TLWP day CAM5','LWP day CAM5'}
            line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='d'; istyle=istyle+1;
        case {'TLWP night CAM5','LWP night CAM5'}
            line_pattern(istyle).p= '--';  line_colour(istyle).c=[0 0 1]; marker_style(istyle).m='d'; istyle=istyle+1;
        case {'TLWP day CAMCLUBBv2','LWP day CAMCLUBBv2'}
            line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 1 0]; marker_style(istyle).m='d'; istyle=istyle+1;
        case {'TLWP night CAMCLUBBv2','LWP night CAMCLUBBv2'}
            line_pattern(istyle).p= '--';  line_colour(istyle).c=[0 1 0]; marker_style(istyle).m='d'; istyle=istyle+1;
    end
end