function [boot_out,boot_out_std]=bootstrap_array(dat,Nboot,prctiles)

%bootstrapping error analysis of each lat/lon cell - sampling n days from
%all that are available.
%Use like this for e.g.:- (for Nd)
%boot_out=bootstrap_array(N_time3,1000,[2.5 5 7.5 10 20 30 50 70 80 90 92.5
%95 97.5]);
%or maybe use dat_modis as from plot_global_maps in order to use data with
%certain data screened out ==>
% boot_out=bootstrap_array(dat_modis,1000,[2.5 5 7.5 10 20 30 50 70 80 90 92.5 95 97.5]);


boot_out = NaN*ones(size(dat,1),size(dat,2),length(prctiles));
boot_out_std = NaN*ones(size(dat,1),size(dat,2));

figure
for ilat=1:size(dat,1)
    fprintf(1,'\n%d',ilat);
    for ilon=1:size(dat,2)

        X2=dat(ilat,ilon,:);

        %Nd_data = squeeze(N_time3(ilat,ilon,:));
        Nd_data = squeeze(X2);
        Nd_data(isnan(Nd_data))='';

        %bootstat = bootstrp(10^4, @(x) [mean(x) std(x)], Nd_data);

        bootstat = bootstrp(Nboot, @(x) [mean(x)], Nd_data);

        pdf_mean = ksdensity(bootstat(:,1)); %pdf of the means of the samples
        %pdf_std = ksdensity(bootstat(:,2)); %pdf of the means of the samples
        
        plot(pdf_mean); hold on


        % prctile(bootstat(:,1),[5 95])
        % 100*diff(prctile(bootstat(:,1),[5 95])) / 2 / mean(bootstat(:,1))
        % prctile(bootstat(:,2),[5 95])
        % mean(bootstat(:,2))
        % 100*mean(bootstat(:,2)) / mean(bootstat(:,1))

        boot_out(ilat,ilon,:) = prctile(bootstat,prctiles);
        
        boot_out_std(ilat,ilon) = std(bootstat);

    end
end




% help bootstrp

% BOOTSTRP Bootstrap statistics.
%     BOOTSTAT = BOOTSTRP(NBOOT,BOOTFUN,D1,...) draws NBOOT bootstrap data
%     samples, computes statistics on each sample using the function BOOTFUN,
%     and returns the results in the matrix BOOTSTATS.  NBOOT must be a
%     positive integer.  BOOTFUN is a function handle specified with @.
%     Each row of BOOTSTAT contains the results of applying BOOTFUN to one
%     bootstrap sample.  If BOOTFUN returns a matrix or array, then this
%     output is converted to a row vector for storage in BOOTSTAT.
%  
%     The third and later input arguments (D1,...) are data (scalars,
%     column vectors, or matrices) that are used to create inputs to BOOTFUN.
%     BOOTSTRP creates each bootstrap sample by sampling with replacement
%     from the rows of the non-scalar data arguments (these must have the
%     same number of rows).  Scalar data are passed to BOOTFUN unchanged.
%  
%     [BOOTSTAT,BOOTSAM] = BOOTSTRP(...) returns BOOTSAM, a matrix of indices
%     into the rows of the extra arguments.  To get the output samples BOOTSAM
%     without applying a function, set BOOTFUN to empty ([]).
%  
%     Examples:
%  
%     Compute a sample of 100 bootstrapped means of random samples taken from
%     the vector Y, and plot an estimate of the density of these bootstrapped
%     means:
%        y = exprnd(5,100,1);
%        m = bootstrp(100, @mean, y);
%        [fi,xi] = ksdensity(m);
%        plot(xi,fi);
%  
%     Compute a sample of 100 bootstrapped means and standard deviations of
%     random samples taken from the vector Y, and plot the bootstrap estimate
%     pairs:
%        y = exprnd(5,100,1);
%        stats = bootstrp(100, @(x) [mean(x) std(x)], y);
%        plot(stats(:,1),stats(:,2),'o')
%  
%     Estimate the standard errors for a coefficient vector in a linear
%     regression by bootstrapping residuals:
%        load hald ingredients heat
%        x = [ones(size(heat)), ingredients];
%        y = heat;
%        b = regress(y,x);
%        yfit = x*b;
%        resid = y - yfit;
%        se = std(bootstrp(1000, @(bootr) regress(yfit+bootr,x), resid));
%  
%     Bootstrap a correlation coefficient standard error:
%        load lawdata gpa lsat
%        se = std(bootstrp(1000,@corr,gpa,lsat));
%  
%     See also random, randsample, hist, ksdensity.
% 
%     Reference page in Help browser
%        doc bootstrp