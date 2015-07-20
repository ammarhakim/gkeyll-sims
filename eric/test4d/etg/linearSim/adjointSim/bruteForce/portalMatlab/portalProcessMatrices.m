% Use to run tests on the cluster
% load data2.mat
function portalProcessMatrices
  parfor i = 1:16
    fprintf('%d\n',i)
  end
end