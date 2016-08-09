clear all; close all;clc

F = importdata('F.dat');
V = importdata('VVEC.dat');

Ftilde2 = importdata('FTILDE-EP2.dat');
JAC2    = importdata('JAC-EP2.dat');

Ftilde4 = importdata('FTILDE-EP4.dat');
JAC4    = importdata('JAC-EP4.dat');

Ftilde6 = importdata('FTILDE-EP6.dat');
JAC6    = importdata('JAC-EP6.dat');

Ftilde8 = importdata('FTILDE-EP8.dat');
JAC8    = importdata('JAC-EP8.dat');

Ftilde10 = importdata('FTILDE-EP10.dat');
JAC10    = importdata('JAC-EP10.dat');

Ftilde12 = importdata('FTILDE-EP12.dat');
JAC12    = importdata('JAC-EP12.dat');

epsilon  = [1e-2;1e-4;1e-6;1e-8;1e-10;1e-12];

norms = [norm(abs(Ftilde2'  -  F' - epsilon(1)*JAC2*V),2);...
         norm(abs(Ftilde4'  -  F' - epsilon(2)*JAC4*V),2);...
         norm(abs(Ftilde6'  -  F' - epsilon(3)*JAC6*V),2);...
         norm(abs(Ftilde8'  -  F' - epsilon(4)*JAC8*V),2);...
         norm(abs(Ftilde10' -  F' - epsilon(5)*JAC10*V),2);...
         norm(abs(Ftilde12' -  F' - epsilon(6)*JAC12*V),2)];

figure; loglog(epsilon,norms,'ko-');ylabel('log(||F(x+\epsilon v)-F(x)-\epsilon J v||)');
xlabel('log (\epsilon)')


