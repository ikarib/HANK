clear; clc; close all
temp = load('TRANS/deltatransvec.txt');
deltatransvec = temp(1,:); worldbond = temp(2,:); rb = temp(3,:);
w = [0,1.7397938975471*(rb(1:199)-0.005)];
for it = 1:199
    w(it+1) = w(it) + 0.1*deltatransvec(it)*(w(it+1)-w(it));
end
max(abs(w-worldbond))