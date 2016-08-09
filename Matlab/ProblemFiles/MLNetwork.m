function F = MLNetwork(t,u,p,W,x,idx)

 %% Rename parameters
  beta = p(16);
  nx   = length(u)/4;
  dx   = x(2)-x(1); 

  %% Splitting
  iV   = idx(:,1); iN = idx(:,2); iC=idx(:,3); iS = idx(:,4);
  v    = u(iV); n = u(iN); c = u(iC);s = u(iS); z = zeros(nx,1);
 
  %% Compute ML neuron component and firing rate
  H = ComputeFiringRate(v,p);
  M = FSMLNeuron(t,u([iV;iN;iC]),p,idx);

  %% Right-hand side
  F   = zeros(4*nx,1);
  F([iV;iN;iC]) = M + [s; z; z];
  F(iS)         =  -beta*s+ (W*H)*dx;

end
