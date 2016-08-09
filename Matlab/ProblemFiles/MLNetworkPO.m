function [G,DGz] = MLNetworkPO(z,p,W,x,idx,zTilde,uTemp,DtuTemp,tTemp)

  %% Split variables
  nU = size(z,1) - 1; iU = 1:nU; iT = nU + 1;
  U = z(iU); T = z(iT);

  %% Problem handles
  problemHandle    = @(t,u) MLNetwork(t,u,p,W,x,idx);
  timeOutputHandle = @(t,u,flag) TimeOutput(t,u,flag);

  %% TimeStep
  tspan   = [0 T];
  options = [];
  [t,UHist]   = ode45(problemHandle,tspan,U,options);

  %% Right-hand side
  G = zeros(size(z));
  G(iU) = U-UHist(end,:)';

  %% Phase condition
  nt = length(tTemp);
  nx = length(x);
  v0 = UHist(:,idx(nx/2,1));
  n0 = UHist(:,idx(nx/2,2));
  c0 = UHist(:,idx(nx/2,3));
  s0 = UHist(:,idx(nx/2,4));
  
  v0 = interp1(t,v0,tTemp);
  n0 = interp1(t,n0,tTemp);
  c0 = interp1(t,c0,tTemp);
  s0 = interp1(t,s0,tTemp);

%   figure; plot(tTemp,v0,tTemp,uTemp(1:nt),'.-')
%   figure; plot(tTemp,n0,tTemp,uTemp(nt   + [1:nt]),'.-')
%   figure; plot(tTemp,c0,tTemp,uTemp(2*nt + [1:nt]),'.-')
%   figure; plot(tTemp,s0,tTemp,uTemp(3*nt + [1:nt]),'.-')

  if any(isnan(v0)) || any(isnan(n0)) || any(isnan(c0)) || any(isnan(s0))
    error('NAN in interpolation');
  end

  G(iT) = DtuTemp'*([v0; n0; c0; s0] - uTemp);

  % PlotHistory(x,t,UHist,p,[],idx,true);

  %% Jacobian-vector action
  if nargout > 1 && ~isempty(zTilde)
    epsi = 1e-4;
    DG  =  MLNetworkPO(z+epsi*zTilde,p,W,x,idx,[],uTemp,DtuTemp,tTemp);
    DGz = (DG - G)/epsi;
  end

 
end
