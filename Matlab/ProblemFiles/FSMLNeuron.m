function M = FSMLNeuron(t,u,p,idx)

  %% Rename parameters
  V1    = p(1) ;
  V2    = p(2) ;
  V3    = p(3) ;
  V4    = p(4) ;
  GCa   = p(5) ;
  GK    = p(6) ;
  GL    = p(7) ;
  VCa   = p(8) ;
  VK    = p(9) ;
  VL    = p(10);
  C     = p(11);
  I0    = p(12);
  IApp  = p(13);
  k     = p(14);
  VStar = p(15);
  rho   = p(24);
  GKCa  = p(25);
  q     = p(26);
  epsilon = p(27);
  mu    = p(28);
  
  nx = size(u,1)/3;

  %% Splitting
  iV   = idx(:,1)'; iN = idx(:,2); iCa=idx(:,3);
  v = u(iV); n = u(iN); Ca= u(iCa);
  
  %% Auxiliary functions
  z      = @(Ca) (Ca)./(Ca + 10);  
  mInf   = @(v) ( 1 + tanh( (v-V1)/V2 ) )/2;
  nInf   = @(v) ( 1 + tanh( (v-V3)/V4 ) )/2;
  r      = @(v) rho*cosh( (v-V3)/(2*V4) );
  ICa    = @(v) GCa*mInf(v).*(v-VCa);
  ITotal = @(v,n,Ca) -(GK*n + GKCa*z(Ca)).*(v - VK) -GL*( v - VL)...
                     - ICa(v) + IApp;

  %% Right-hand side
  M = zeros(size(u));
  M(iV) = ITotal(v,n,Ca)/C;
  M(iN) = r(v) .* ( nInf(v) - n);
  M(iCa)= epsilon*(-mu*ICa(v) - Ca);

    
end
