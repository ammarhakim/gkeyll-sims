function processMatrices_3
  load innerProductMatrix3.mtx
  
  % Set value at third column, first row to 0.0, which is usually the number
  % of nonzero entries in the sparse matrix but matlab will interpret it as
  % an entry if not set to zero
  innerProductMatrix3(1,3) = 0.0;
  sparseM = spconvert(innerProductMatrix3);
  spy(sparseM)
  L = chol(sparseM,'lower'); % Returns lower triangular matrix satisfying
  % L*L' = sparseM
  F = L.';
  
  clear innerProductMatrix3
  clear L
  clear sparseM
  
  load linearOperatorMatrix3.mtx
  linearOperatorMatrix3(1,3) = 0.0;
  sparseLO = spconvert(linearOperatorMatrix3);
  clear linearOperatorMatrix3
  % spy(sparseLO)
  [V,D] = eig(full(sparseLO));
  clear sparseLO
  
  B = F*V;
  
  % Save workspace
  save('matlab_3.mat','-v7.3')
end
