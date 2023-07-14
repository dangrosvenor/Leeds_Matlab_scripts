ichoose_styles=1;

for idat2=1:length(labs)
    switch labs(idat2).l
        case {'TLWP day AMSRE','LWP day AMSRE'}
            line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
            labs(idat2).l = 'AMSRE day';
        case {'TLWP night AMSRE','LWP night AMSRE'}
            line_pattern(istyle).p= '--';  line_colour(istyle).c=[1 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
            labs(idat2).l = 'AMSRE night';            
        case {'TLWP day CAM5','LWP day CAM5'}
            line_pattern(istyle).p= '-';  line_colour(istyle).c=[1 0.65 0]; marker_style(istyle).m='d'; istyle=istyle+1;
            labs(idat2).l = 'CAM5 day';
        case {'TLWP night CAM5','LWP night CAM5'}
            line_pattern(istyle).p= '--';  line_colour(istyle).c=[1 0.65 0]; marker_style(istyle).m='d'; istyle=istyle+1;
            labs(idat2).l = 'CAM5 night';            
        case {'TLWP day CAMCLUBBv2','LWP day CAMCLUBBv2'}
            line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 1 0]; marker_style(istyle).m='d'; istyle=istyle+1;
            labs(idat2).l = 'CAMCLUBBv2 day';            
        case {'TLWP night CAMCLUBBv2','LWP night CAMCLUBBv2'}
            line_pattern(istyle).p= '--';  line_colour(istyle).c=[0 1 0]; marker_style(istyle).m='d'; istyle=istyle+1;
            labs(idat2).l = 'CAMCLUBBv2 night';             
        case {'LWP day AM3'}
            line_pattern(istyle).p= '-';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='d'; istyle=istyle+1;
            labs(idat2).l = 'AM3 day';             
        case {'LWP night AM3'}
            line_pattern(istyle).p= '--';  line_colour(istyle).c=[0 0 0]; marker_style(istyle).m='d'; istyle=istyle+1; 
            labs(idat2).l = 'AM3 night';                         
        case {'LWP day AM3 CLUBBv1 2deg'}
            line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.4 0.4 0.4]; marker_style(istyle).m='d'; istyle=istyle+1;
            labs(idat2).l = 'AM3CLUBBv1 day';                         
        case {'LWP night AM3 CLUBBv1 2deg'}
            line_pattern(istyle).p= '--';  line_colour(istyle).c=[0.4 0.4 0.4]; marker_style(istyle).m='d'; istyle=istyle+1;            
            labs(idat2).l = 'AM3CLUBBv1 night';                
        case {'LWP day AM3 CLUBBv2 COSP 200km'}
            line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.4 0.4 1]; marker_style(istyle).m='d'; istyle=istyle+1;
            labs(idat2).l = 'AM3CLUBBv2 day';                
        case {'LWP night AM3 CLUBBv2 COSP 200km'}
            line_pattern(istyle).p= '--';  line_colour(istyle).c=[0.4 0.4 1]; marker_style(istyle).m='d'; istyle=istyle+1;  
            labs(idat2).l = 'AM3CLUBBv2 night';                
        case {'LWP day CAM5 CLUBB COSP'}  %[0 0 1]
            line_pattern(istyle).p= '-';  line_colour(istyle).c=[0.7 0.4 0.7]; marker_style(istyle).m='d'; istyle=istyle+1;
            labs(idat2).l = 'CAM5CLUBBv1 day';                
        case {'LWP night CAM5 CLUBB COSP'}
            line_pattern(istyle).p= '--';  line_colour(istyle).c=[0.7 0.4 0.7]; marker_style(istyle).m='d'; istyle=istyle+1;     
            labs(idat2).l = 'CAM5CLUBBv1 night';                            
        
            
    end
end