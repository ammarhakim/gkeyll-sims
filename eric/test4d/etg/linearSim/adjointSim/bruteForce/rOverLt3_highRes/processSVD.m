function processMatrices_3
  data = load('matlab_3.mat','F','V','D','B');
  F = data.F;
  V = data.V;
  D = data.D;
  B = data.B;
  clear data

  gammaMax = max(real(diag(D)));
  tPoints = linspace(0, 40e-6, 10);
  gPoints = zeros(1,length(tPoints));

  parfor i = 1:length(tPoints)
    t = tPoints(i);
     [U_svd,S_svd,V_svd] = svdsecon(F*V*diag(exp(diag(D*t)))/B,1);
     gPoints(i) = S_svd(1,1)^2;
    fprintf('G = %d, t = %d\n', gPoints(i),tPoints(i));
  end
  
  % Write to file
  plotData(:,1) = tPoints';
  plotData(:,2) = gPoints';
  save('data3.mat','plotData','gammaMax')
end
