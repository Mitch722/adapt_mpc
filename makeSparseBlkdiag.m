function [blkA_sparse, blkB_sparse, blkC_sparse, kalm_gain] = makeSparseBlkdiag(hyper_params, no_states, Ts)


[no_models, ~] = size(hyper_params);

blkA = zeros(no_models*no_states);
blkBcell = cell(no_models, 1);
blkC = zeros(no_models*2, no_models*no_states);
blkL = zeros(4*no_models, 2*no_models);

for i = 1 : no_models
    
    j = i-1;
    
    a11 = hyper_params(i, 1);
    a22 = hyper_params(i, 2);
    
    b11 = hyper_params(i, 3);
    b22 = hyper_params(i, 4);
    
    obvs_poles = [0.01, 0.1, 0.02, 0.02]';
    
    A = makesysd_a(a11, a22, b11, b22, Ts);
    
    L = place(A.A', A.C', obvs_poles);
    L = L';
    
    blkA(j*no_states+1 : (j+1)*no_states, j*no_states+1 : (j+1)*no_states) = A.A;
    blkBcell{i, 1} = A.B;
    blkC(j*2+1 : (j+1)*2, j*no_states+1 : (j+1)*no_states) = A.C;
    blkL(j*4+1 : (j+1)*4, j*2+1 : (j+1)*2) = L;
    
    
end


blkA_sparse = sparse(blkA);
blkB_sparse = sparse(cell2mat(blkBcell));
blkC_sparse = sparse(blkC);
kalm_gain = sparse(blkL);

