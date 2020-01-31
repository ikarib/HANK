%% This function is doing one of these
% parfor iy = 1:ngpy
%     x(:,:,iy) = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)'\x(:,:,iy); % or
%     x(:,iy) = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)'\x(:,iy); % or
%     x(:,iy) = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)\x(:,iy); % or
%     x(:,iy) = spdiags(B(:,:,iy,it),[0 -1 1 -ngpa ngpa],nab,nab)\x(:,iy);
% end
% syntax:
%   spsolve(transb,ngpa,B,x)
% where transb is 'T' for transpose or 'N' for no transpose

function x = spsolve(transb,ngpa,B,x)
assert(size(B,2)==5)
[nab,~,ngpy] = size(B);
assert(size(x,1)==nab)
use_umfpack = true;
for iy = 1:ngpy %par
    if size(x,2)==ngpy
        if transb=='T'
            if use_umfpack
                x(:,iy)   = umfpack(spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)','\',x(:,iy));
            else
                x(:,iy)   = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)'\x(:,iy);
            end
        elseif transb=='N'
            if use_umfpack
                x(:,iy)   = umfpack(spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab),'\',x(:,iy));
            else
                x(:,iy)   = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)\x(:,iy);
            end
        end
    elseif size(x,3)==ngpy
        if transb=='T'
%             if use_umfpack
%                 x(:,:,iy) = umfpack(spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)','\',x(:,:,iy));
%             else
                x(:,:,iy) = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)'\x(:,:,iy);
%             end
        elseif transb=='N'
            if use_umfpack
                x(:,:,iy) = umfpack(spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab),'\',x(:,:,iy));
            else
                x(:,:,iy) = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)\x(:,:,iy);
            end
        end
    end
end