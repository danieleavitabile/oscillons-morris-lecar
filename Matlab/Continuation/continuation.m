% Cleaning 
clear all, close all, clc;

% Load problem files and solution
addpath('~/Dropbox/SecantContinuation/0.2/src');
addpath('../ProblemFiles');
addpath('../SolutionsDatabase');

%% Define System 
nx  = 460; x   = linspace(0,60,nx)'; hx  = abs(x(1)-x(2));
iV  = (1:nx)'; iN = nx+iV; iC = nx+iN; iS = nx+iC;
idx = [iV iN iC iS];

%% Initial guess
sol = load('initialOscillon');
z0 = sol.u;
p0 = sol.p;

%% Connectivity Matrix
[W,~] = SynapticKernel(p0,x);

%% Determine template for phase condition
problemHandle    = @(t,u) MLNetwork(t,u,p,W,x,idx);
timeOutputHandle = @(t,u,flag) TimeOutput(t,u,flag);

%% TimeStep
nT  = 100; tspan   = linspace(0,z0(end),nT+1); ht = tspan(2) - tspan(1); iT = 1:nT;
[tTemp,UHist]   = ode45(@(t,u) MLNetwork(t,u,p0,W,x,idx),tspan,z0(1:end-1));

e = ones(nT,1);
Dt = spdiags([-e e],[-1 1],nT,nT); Dt(1,nT) = -1; Dt(nT,1) = 1; Dt = Dt/ht;
uTemp   = [UHist(iT,iV(nx/2)); UHist(iT,iN(nx/2)); UHist(iT,iC(nx/2)); UHist(iT,iS(nx/2)) ];
DtuTemp = [Dt*UHist(iT,iV(nx/2)); Dt*UHist(iT,iN(nx/2)); Dt*UHist(iT,iC(nx/2)); Dt*UHist(iT,iS(nx/2)) ];
% plot(tTemp(iT),uTemp(iT),tTemp(iT),DtuTemp(iT));

prob    = @(z,p) MLNetworkPO(z,p,W,x,idx,uTemp,DtuTemp,tTemp(iT));
plotSol = @(z,p,parent) PlotSolution(x,z(1:end-1),p,parent,idx,false);
F = prob(z0,p0);


% plotSol(F,p0,[]);
% 
% % 
% % NeuralFieldForcedPO(u,p,wHat,x,Lx,idx,[]);
% % jac      = @(u,p,v) NeuralFieldForcedPOJacobianAction(u,p,wHat,x,Lx,idx,v);
% % plotSol  = @(u,p,parent) PlotSolutionPO(u,p,parent,wHat,x,Lx,idx,true); % PlotSolution(x,u,p,parent,idx,false);
% % solMeas  = @(step,u,p) SolutionMeasures(step,u,p,x,idx,'periodic');
% % compSpec = [];%@(u,p) ComputeSpectrum(u,p,wHat,x,Lx,idx);
% % plotSpec = [];%@(d,p,parent) PlotSpectrum(d,p,parent);
% % 
% % %% Assign problem 
% % stepPars.iContPar                         = 3;
% % stepPars.pMin                             = -60;
% % stepPars.pMax                             = 60;
% % stepPars.s0                               = -0.8;
% % stepPars.sMin                             =  0.008;
% % stepPars.sMax                             =  1.0;
% % stepPars.maxSteps                         = 20000;
% % stepPars.nPrint                           = 1;
% % stepPars.nSaveSol                         = 1;
% % stepPars.finDiffEps                       = 1e-5;
% % stepPars.optNonlinIter                    = 5;
% % stepPars.NewtonGMRESOptions.nonlinTol     = 1e-7;
% % stepPars.NewtonGMRESOptions.nonlinMaxIter = 7;
% % stepPars.NewtonGMRESOptions.linTol        = 1e-3;
% % stepPars.NewtonGMRESOptions.linrestart    = 1000;
% % stepPars.NewtonGMRESOptions.linmaxit      = 50;
% % stepPars.NewtonGMRESOptions.damping       = 1.0;
% % stepPars.NewtonGMRESOptions.display       = 1;
% % stepPars.dataFolder                       = 'Data';
% % stepPars.PlotSolution                     = plotSol;
% % stepPars.BranchVariables                  = solMeas;
% % stepPars.PlotBranchVariableId             = 4;
% % stepPars.ComputeEigenvalues               = compSpec;
% % stepPars.PlotSpectrum                     = plotSpec;
% % 
% % %% Run
% % branch = SecantContinuationNewtonGMRES(prob,jac,u0,p0,stepPars);
