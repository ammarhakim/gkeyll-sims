% Use to run tests on the cluster
% load data2.mat
function portalProcessMatrices
  data = load('matlab_2.mat','F','V','D','B','gammaMax');
  F = data.F;
  V = data.V;
  D = data.D;
  B = data.B;
  gammaMax = data.gammaMax;
  clear data
  
  tPoints = linspace(0, 40e-6, 16);
  gPoints = zeros(1,length(tPoints));

  parfor i = 1:length(tPoints)
    t = tPoints(i);
    [U_svd,S_svd,V_svd] = svdsecon(F*V*diag(exp(diag(D*t)))/B,1);
    gPoints(i) = S_svd(1,1)^2;
    fprintf('G = %d, t = %d\n', gPoints2(i),tPoints2(i));
  end
  
  % Write to file
  plotData(:,1) = tPoints';
  plotData(:,2) = gPoints';
  save('data2.mat','plotData','gammaMax')
end