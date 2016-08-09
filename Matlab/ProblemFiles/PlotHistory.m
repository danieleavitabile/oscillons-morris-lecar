function plotHandle = PlotHistory(x,t,U,p,parentHandle,idx,plotTraces)

   %% Position and eventually grab figure
   if isempty(parentHandle)
     scrsz = get(0,'ScreenSize');
     plotHandle = figure('Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/4 scrsz(4)/4]);
     parentHandle = plotHandle;
    else
      plotHandle = parentHandle;
   end
   figure(parentHandle);

   %% Extract number of components
   numComp = size(idx,2);

   %% Grid
   [T,X] = meshgrid(t,x);
   nx = length(x);

   %% Plots
   if ~plotTraces
     for k = 1:numComp
       subplot(1,numComp,k)
       pcolor(X,T,U(:,idx(:,k))'); shading interp; view([0 90]);
       title(['U_' num2str(k)]);
       xlabel('x'); ylabel('t');
     end
   else

     for k = 1:numComp
       subplot(2,numComp,k)
       pcolor(X,T,U(:,idx(:,k))'); shading interp; view([0 90]);
       title(['U_' num2str(k)]);
       xlabel('x'); ylabel('t');
     end

     subplot(2,numComp,[numComp+1 : 2*numComp]);% hold on;
     u10 = U(:,idx(nx/2,1)); u20 = U(:,idx(nx/2,2));
     u30 = U(:,idx(nx/2,3)); u40 = U(:,idx(nx/2,4));
     plot(t,[u10 u20 u30 ,u40],'.-');
     xlim([min(t) max(t)]);
     xlabel('t'); 
     legend('R(0,t)','V(0,t)','R(L,t)','V(L,t)')

   end

   %% Save
   print -dtiff history.tiff

end

