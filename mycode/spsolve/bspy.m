clear;load spsolve7
ngpb=50;
%%
for iy=1:ngpy
    subplot(1,2,1)
    b=B(:,2,iy,1)+B(:,3,iy,1);
    surf(reshape(b,ngpa,ngpb))

    subplot(1,2,2)
    b=B(:,4,iy,1)+B(:,5,iy,1);
    surf(reshape(b,ngpa,ngpb))
    pause(1)
end