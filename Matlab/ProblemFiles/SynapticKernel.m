function [W,wHandle] = SynapticKernel(p,x)

  %% Parameters
  A1    = p(18);
  B1    = p(20);
  L     = p(23);
  
  %% Function handles
  %   wf = @(z,A,B) A/sqrt(pi*B) * exp(-(z.^2)/ B);   %Gaussian
  wf = @(z,A,B) A * exp(-B*abs(z));
  wh = @(z) wf(z,A1,B1);                              %Exponential 
  wHandle  = @(z) wh(x-L/2);  
 
  %% Build coupling matrix
  N = length(x);
  j = 1:N;
  i = N/2;
  W = zeros(N,N);

  W(N/2,:)=N*wh(2*L*abs(i-j)/N);
  for i = N/2+1:N
    W(i,:) = circshift( W(i-1,:)',  1 )';
  end
  for i = N/2-1:-1:1
    W(i,:) = circshift( W(i+1,:)', -1 )';
  end
   W = W/N;
  

end
