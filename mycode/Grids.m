global agrid bgrid ygrid dagrid dbgrid ymarkovdiag ymarkovoff

% productivity process
if ngpy==1
	logygrid = 0;
	ygrid = 1;
	ymarkov = 0;
	ytrans = 1;
	ydist = 1;
elseif ReadEarningsProcess
    fid1 = fopen(fullfile(EarningsProcessDir,'ygrid_combined.txt'),'r'); logygrid = fscanf(fid1,'%f'); fclose(fid1);
	fid1 = fopen(fullfile(EarningsProcessDir,'ydist_combined.txt'),'r'); ydist = fscanf(fid1,'%f'); fclose(fid1);
    fid1 = fopen(fullfile(EarningsProcessDir,'ymarkov_combined.txt'),'r'); ymarkov = fscanf(fid1,'%f'); fclose(fid1);
	if AdjustProdGridFrisch; logygrid = logygrid/(1+adjfricshgridfrac*frisch); end
	ygrid = exp(logygrid);
    ymarkov = reshape(ymarkov,ngpy,ngpy)'; % since fortran reads in column major order
    ymarkov = ymarkov-diag(sum(ymarkov,2)); % fix up rounding in markov matrix
    ydist = ydist/sum(ydist); % fix up rounding in ergodic distribution
elseif TwoPointWageProcess
	ygrid = [0.8 1.2]';
	logygrid = log(ygrid);
	ytrans = [1-0.06667 0.06667; 0.06667 1-0.06667];
	ydist = [0.5 0.5]';
	ymarkov = ytrans-eye(ngpy); % assumes ytrans is quarterly	
end

ymarkovdiag = diag(ymarkov);
ymarkovoff = ymarkov-diag(ymarkovdiag);
ymarkovdiag = shiftdim(ymarkovdiag,-2);

% adjust mean productivity
lmean = ydist'*ygrid;
ygrid = shiftdim(meanlabeff*ygrid/lmean,-2);
ydist = shiftdim(ydist,-2);

agrid = PowerSpacedGrid(ngpa,agridparam,0,amax)';

% with low gridparam points get bunched close to zero, so evenly space first 8 points;
if ngpa>10
	agrid(1:9) = (0:8)*agrid(10)/9;
end
dagrid = diff(agrid);

% liquid grid
if Borrowing
	nbl = -lumptransfer/(rborr+PerfectAnnuityMarkets*deathrate);
 	abl = max(nbl+cmin,blim);
    if Display>1
        fprintf('natural borrowing limit = %.15g\n',nbl);
        fprintf('actual borrowing limit = %.15g\n',abl);
    end
    bgrid = PowerSpacedGrid(ngpbNEG/2+1,bgridparamNEG,abl,abl/2);
	bgrid = [bgrid abl-fliplr(bgrid(2:end-1)) PowerSpacedGrid(ngpbPOS,bgridparam,0,bmax)];
else
	bgrid = PowerSpacedGrid(ngpb,bgridparam,0,bmax);
end
dbgrid = diff(bgrid);

diff2 = @(X) (X([2:end end])-X([1 1:end-1]))/2;
adelta = diff2(agrid);
bdelta = diff2(bgrid);
abdelta = adelta*bdelta;
ABdelta = abdelta(:)'./abdelta(:);
ABdelta = [[diag(ABdelta,-1);0] ...
           [0;diag(ABdelta,1)] ...
           [diag(ABdelta,-ngpa);zeros(ngpa,1)] ...
           [zeros(ngpa,1);diag(ABdelta,ngpa)]];
