% !nvcc simulate.cu -lknitro -DDEBUG
% !srun --gres=gpu:1 ./a.out
clear;clc;close all
load obj.txt
[objmin,i_min]=min(obj);
lambda=0;sigma=1.5; delta=0:.001:1;
siz=[length(delta),1,length(sigma),1,length(lambda),1];
[Dy,Dx,Sy,Sx,Ly,Lx]=ind2sub(siz,i_min);
obj=reshape(obj,siz);
fprintf('lambda\t%.15g\t%.15g\n',lambda(Lx),lambda(Ly));
fprintf('sigma\t%.15g\t%.15g\n',sigma(Sx),sigma(Sy));
fprintf('delta\t%.15g\t%.15g\n',delta(Dx),delta(Dy));
fprintf('objfun\t%.15g\t%.15g\n',obj(Dy,Dx,Sy,Sx,Ly,Lx),objmin)
subplot(2,3,1); plot(lambda,squeeze(obj(Dy,Dx,Sy,Sx,Ly,:))); xlabel('\lambda_1')
subplot(2,3,2); plot(sigma,squeeze(obj(Dy,Dx,Sy,:,Ly,Lx))); xlabel('\sigma_1')
subplot(2,3,3); plot(delta,squeeze(obj(Dy,:,Sy,Sx,Ly,Lx))); xlabel('\delta_1')
subplot(2,3,4); plot(lambda,squeeze(obj(Dy,Dx,Sy,Sx,:,Lx))); xlabel('\lambda_2')
subplot(2,3,5); plot(sigma,squeeze(obj(Dy,Dx,:,Sx,Ly,Lx))); xlabel('\sigma_2')
subplot(2,3,6); plot(delta,squeeze(obj(:,Dx,Sy,Sx,Ly,Lx))); xlabel('\delta_2')