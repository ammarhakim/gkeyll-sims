% Use to run tests on the cluster
% load data2.mat
matlabpool open local 8

parfor i = 1:16
  fprintf('%d\n',i)
end

matlabpool close

exit