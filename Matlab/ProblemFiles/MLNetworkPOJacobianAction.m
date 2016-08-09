function DGz = MLNetworkPOJacobianAction(z,p,W,x,idx,zTilde,uTemp,DtuTemp,tTemp);
  [~,DGz] = MLNetworkPO(z,p,W,x,idx,zTilde,uTemp,DtuTemp,tTemp);
end
