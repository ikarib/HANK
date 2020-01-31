clear
N=20000;
A=rand(N);
B=rand(N);
disp('done')
%%
tic;C=max(A,B);toc
tic;
for i=1:100
    max_in_place(A,B);
end
toc
assert(all(all(A==C)))