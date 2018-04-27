% spsolve(x,B(:,:,:,it))
% spsolve(x,B)
% spsolve(lvec,B)
% spsolve(rhs,B)
% spsolve(x,B)
function x = spsolve(B,x)

if size(x,3)>1
    for iy=1:ngpy %par
        x(:,:,iy) = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)'\x(:,:,iy);
    end
else
    for iy=1:ngpy %par
        x(:,iy)   = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)\x(:,iy);
    end
end