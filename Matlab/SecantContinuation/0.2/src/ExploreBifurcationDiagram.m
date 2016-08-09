function plotHandle = ExploreBifurcationDiagram(branchFile,idVar)

  if isempty(branchFile)
    [fileName,pathName] = uigetfile('*.mat','Select the branch file');
    branchFile = [pathName fileName];
  end

  branchData = load(branchFile);
  branch = branchData.branch;
  plotHandle = figure;

  hold on;
  plot(branch(:,3),branch(:,idVar),'b-')

  iUnstab = find(branch(:,2) > 0);
  iStab   = find(branch(:,2) == 0); 
  iNoStab = find(branch(:,2) == -1); 
  if ~isempty(iStab)
    plot(branch(iStab,3),branch(iStab,idVar),'b.');
  end
  if ~isempty(iUnstab)
    plot(branch(iUnstab,3),branch(iUnstab,idVar),'r.');
  end
  if ~isempty(iNoStab)
    plot(branch(iNoStab,3),branch(iNoStab,idVar),'k.');
  end

  numPoints = size(branch,1);
  for i = 1:40:numPoints
    text(branch(i,3),branch(i,idVar),num2str(i-1));
  end
  hold off;

end
