clear all; close all;clc;

%% Load from Matlab 
addpath('../ProblemFiles/');
addpath('../SolutionsDatabase/');

sol    = load('MatlabOscillon-00.mat');
p0     = sol.p;
u0     = sol.u;
Period = sol.T;

%% Load from AUTO  solution
% Z      = importdata('AUTOscillon-00.dat');
% u0     = Z(end,:);
% Period = Z(end,1);

%% Define System 
nx  = 460; x   = linspace(0,60,nx)'; hx  = abs(x(1)-x(2));

iV  = (1:nx)'; iN = nx+iV; iC = nx+iN; iS = nx+iC;
idx = [iV iN iC iS];

%% Connectivity Matrix
[W,~] = SynapticKernel(p0,x);

%% Problem handles
problemHandle    = @(t,u) MLNetwork(t,u,p0,W,x,idx);
timeOutputHandle = @(t,u,flag) TimeOutput(t,u,flag);
plotSol          = @(u,p,parent) PlotSolution(x,u,p,parent,idx,false)

%% TimeStep
tspan   = [0 200];
options = odeset('OutputFcn',timeOutputHandle);
[t,U]   = ode45(problemHandle,tspan,u0,options);

plotSol(U(1,:)',p0,[]);
plotSol(U(end,:)',p0,[]);
PlotHistory(x,t,U,p0,[],idx,true);

% Calculate period
[zmax,Ix] = findpeaks(U(:,iV(nx/2)));
Periods   = diff(t(Ix(2:end-1)))
Period    = mean(Periods)

tspan = [0:0.01:12*Period];
[t,U]   = ode45(problemHandle,tspan,U(end,:)',options);

plotSol(U(1,:)',p0,[]);
plotSol(U(end,:)',p0,[]);
% PlotHistory(x,t,U,p0,[],idx,true);

% Calculate period
[zmax,Ix] = findpeaks(U(:,iV(nx/2)));
Periods   = diff(t(Ix(2:end-1)))
Period    = mean(Periods)

% Save final state oscillon for Matlab 
u       = [U(end,:)'; Period];
p       = p0;
save('initialOscillon.mat','u','p');


% %% Manage Output
% 
% % Calculate period
% [zmax,Ix] = findpeaks(U(:,round(iV(nx/2))));
% Periods   = diff(t(Ix(2:end-1)));
% Period    = mean(Periods);
% 
% % Save single oscillon history for AUTO
% uHist   = U(Ix(2):Ix(3),:);
% tHist   = t(Ix(2):Ix(3)) - t(Ix(2));
% Z       = [tHist,uHist];
% save('AUTOscillon-00.dat','Z','-ascii','-double')
% 
% % Save final state oscillon for Matlab 
% u       = U(end,:);
% p       = p0;
% T       = Period;
% save('MatlabOscillon-01.mat','u','p','T')
% 
% %% Plot History
% V = U(:,iV); N = U(:,iN); C = U(:,iC); S = U(:,iS);
% [X,T] = meshgrid(x,t);
% figure;
% subplot(2,3,[1 4]); pcolor(X,T,V); shading interp; 
% subplot(2,3,[2 5]); pcolor(X,T,C); shading interp; 
% subplot(2,3,[3 6]); pcolor(X,T,S); shading interp; 
% 
% %% Plot final profile
% figure;
% subplot(4,1,1);
% plot(x,U(end,iV),'-');
% subplot(4,1,2);
% plot(x,U(end,iN),'-');
% subplot(4,1,3);
% plot(x,U(end,iC),'-');
% subplot(4,1,4);
% plot(x,U(end,iS),'-');
% 
% %% Watch Oscillon Movie
% % tl = length(t);
% % for i =1:tl 
% % z = U(i,:); 
% % h = figure(3);drawnow;
% % subplot(4,1,1);
% % title(['t= ',num2str(t(i))]);
% % plot(x,U(i,iV),'-');xlim([0 x(end)]);ylim([-40 40]);drawnow;
% % subplot(4,1,2);
% % title(['t= ',num2str(t(i))]);
% % plot(x,U(i,iN),'-');xlim([0 x(end)]);ylim([0 1]);drawnow;
% % subplot(4,1,3);
% % title(['t= ',num2str(t(i))]);
% % plot(x,U(i,iC),'-');xlim([0 x(end)]);ylim([20 80]);drawnow;
% % subplot(4,1,4); 
% % title(['t= ',num2str(t(i))]);
% % plot(x,U(i,iS),'-');xlim([0 x(end)]);ylim([0 10]);drawnow;
% % pause(0.01);
% % end
% 
% 
% 
% tspan = [0 T];
% [t,U]   = ode45(problemHandle,tspan,u,options);
% 
% %% Plot History
% V = U(:,iV); N = U(:,iN); C = U(:,iC); S = U(:,iS);
% [X,T] = meshgrid(x,t);
% figure;
% subplot(2,3,[1 4]); pcolor(X,T,V); shading interp; 
% subplot(2,3,[2 5]); pcolor(X,T,C); shading interp; 
% subplot(2,3,[3 6]); pcolor(X,T,S); shading interp; 
% 
