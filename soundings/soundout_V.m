    function []=soundout(p,savedir)

    fid=fopen(savedir,'wt');
	for i=1:size(p,1)
        fprintf(fid,'%g %g %g %g %g %g %g %g\n',999,999,p(i,1),p(i,2),p(i,3),p(i,4),999,p(i,9));
	end
    fclose(fid);