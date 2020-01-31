%% Transition
clear;load spsolve7
x=lgmat1;
tic
for it=1:200
    x=spsolve('N',ngpa,B(:,:,:,it),x);
end
t1=toc
tic
for it=1:200
    for iy = 1:ngpy %par
        lgmat1(:,iy) = spdiags(B(:,:,iy,it),[0 -1 1 -ngpa ngpa],nab,nab)\lgmat1(:,iy);
    end
end
t2=toc
max(max(abs(x-lgmat1)))
fprintf('spsolve took %f seconds and is faster by %f seconds\n',t1,t2-t1)
%% HJBUpdate
clear;load spsolve1 % spsolve6
x=lbvec;
tic
x=spsolve('T',ngpa,B,x);
t1=toc;
tic
for iy = 1:ngpy %par
    lbvec(:,iy) = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)'\lbvec(:,iy);
end
t2=toc;
max(max(abs(x-lbvec)))
fprintf('spsolve took %f seconds and is faster by %f seconds\n',t1,t2-t1)
%% StationaryDistribution
clear;load spsolve2
tic
x=spsolve('N',ngpa,B,lgmat);
t1=toc;
tic
for iy = 1:ngpy %par
    lgmat(:,iy) = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)\lgmat(:,iy);
end
t2=toc;
max(max(abs(x-lgmat)))
fprintf('spsolve took %f seconds and is faster by %f seconds\n',t1,t2-t1)
%% CumulativeConsumption
clear;load spsolve3
tic
x=spsolve('T',ngpa,B,cdcumvec);
t1=toc;
tic
for iy = 1:ngpy %par
    cdcumvec(:,:,iy) = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)'\cdcumvec(:,:,iy);
end
t2=toc;
max(max(max(abs(x-cdcumvec))))
fprintf('spsolve took %f seconds and is faster by %f seconds\n',t1,t2-t1)
%% DiscountedMPC
clear;load spsolve4 % spsolve5,spsolve8
tic
x=spsolve('T',ngpa,B,lvec);
t1=toc;
tic
for iy = 1:ngpy %par
    lvec(:,iy) = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)'\lvec(:,iy);
end
t2=toc;
max(max(abs(x-lvec)))
fprintf('spsolve took %f seconds and is faster by %f seconds\n',t1,t2-t1)