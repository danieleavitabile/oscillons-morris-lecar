
function [S,dS] = ComputeFiringRate(v,p)

   k     = p(14);
   vStar = p(15);

   S = 1 ./ (1 + exp(- k * (v - vStar) ));

   if nargout > 1
      dS = k * S .* (1-S);
   end

end
