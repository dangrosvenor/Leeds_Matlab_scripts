function [me_out, std_out]=std_combine(stds,means,Ns)
%% See std_combine2 now!
%function [me_out, std_out]=std_combine(stds,means,Ns)
%Combine the stds (and means) from several different samples
% stds, means and Ns are vectors of length N of the std devs, means and sample sizes
% for the N populations.
%Uses the formulae from Wikipedia - http://en.wikipedia.org/wiki/Standard_deviation#Combining_standard_deviations

fprintf(1,'\n\n***Superseeded by std_combine2!!! ***\n\n');

Ntot = sum(Ns);
%Ntot = sum(Ns,2);

%Mean
me_out = sum(Ns.*means) ./ Ntot;
%me_out = sum(Ns.*means,2) ./ Ntot;

%Std dev
std_out = sqrt(  sum(Ns.*(stds.^2 + means.^2 - me_out.^2) ) / Ntot  );
%std_out = sqrt(  sum(Ns.*(stds.^2 + means.^2 - me_out.^2) ,2) / Ntot  );

%check - gives the same answer, so just use the simpler formuation above
% A = sum(Ns.*stds.^2)/Ntot;
%
% B=0;
% for i=2:length(Ns)
%     for j=1:i
%         B = B + Ns(i)*Ns(j)*(means(i) - means(j)).^2;
%     end
% end
%
% std_out2 = sqrt( A + B/Ntot.^2 );

%Another check using actual data samples - gives slightly different stds to the true answer (stdA) -
% s1 = 27.86, stdA = 27.83. But probably close enough!
run_test=0;
if run_test==1
    
    test_case=2;

    switch test_case
        case 1

            clear a means stds Ns
            a{1}=magic(5);
            a{2}=magic(10)
            a{3}=magic(7);

            A=[];
            for i=1:3
                means(i) = mean(a{i}(:));
                stds(i) = std(a{i}(:));
                Ns(i) = length(a{i}(:));
                A = cat(1,A,a{i}(:));
            end

            [me,s1] = std_combine(stds,means,Ns)
            meA = mean(A)
            stdA = std(A)

        case 2
            d3 = rand([360,180,12]);
            [me3,N3,std3] = meanNoNan(d3(:),1);
            
            [me2,N2,std2] = meanNoNan(d3,1);
            me1 = meanNoNan(
            


    end

end





