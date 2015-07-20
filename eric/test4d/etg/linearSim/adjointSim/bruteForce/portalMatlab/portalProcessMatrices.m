% Use to run tests on the cluster
% load data2.mat
function portalProcessMatrices
  matlabpool open local 8
  
  disp('hello')
  for i = 1:16
    fprintf('%d\n',i)
  end
  
  matlabpool close
end