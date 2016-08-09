function F = SolutionMeasures(step,u,p)

%   %% Rename parameters
%   ic = length(u);

  %% Assign branch variables
%   F = norm(u(1:end-1),2);
   F = u(end);
  
end
