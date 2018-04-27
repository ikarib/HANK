if Display; fprintf('Computing cumulative consumption for MPC\n'); end

B = [1+deltacumcon*(sum(AU,2)-ymarkovdiag) -deltacumcon*AU];

nsteps = round(4/deltacumcon);
cdvec = reshape([c d],nab,2,ngpy);
cdcumvec = zeros(nab,2,ngpy);
tic
for it = 1:nsteps
    cdcumvec = cdcumvec + deltacumcon*(cdvec + reshape(reshape(cdcumvec,2*nab,ngpy)*ymarkovoff',nab,2,ngpy));
    for iy=1:ngpy %par
        cdcumvec(:,:,iy) = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)'\cdcumvec(:,:,iy);
    end
    if it==round(1/deltacumcon)
        % one quarter cumulative expected consumption and deposits
        ccum1 = reshape(cdcumvec(:,1,:),ngpa,ngpb,ngpy); % consumption
        dcum1 = reshape(cdcumvec(:,2,:),ngpa,ngpb,ngpy); % deposits
    end
    if it==round(2/deltacumcon)
        % two quarter cumulative expected consumption
        ccum2 = reshape(cdcumvec(:,1,:),ngpa,ngpb,ngpy); % consumption
        dcum2 = reshape(cdcumvec(:,2,:),ngpa,ngpb,ngpy); % deposits
    end
end
toc
% four quarter cumulative expected consumption
ccum4 = reshape(cdcumvec(:,1,:),ngpa,ngpb,ngpy); % consumption
dcum4 = reshape(cdcumvec(:,2,:),ngpa,ngpb,ngpy); % deposits

cumINITSS = struct('ccum1',ccum1,'ccum2',ccum2,'ccum4',ccum4,...
                   'dcum1',dcum1,'dcum2',dcum2,'dcum4',dcum4);

clear iy nsteps cdvec cdcumvec