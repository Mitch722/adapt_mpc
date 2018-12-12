function blkA_sparse = makeSparseBlkdiag(hyper_params, no_states, Ts)


[no_models, ~] = size(hyper_params);

blkA = zeros(no_models*no_states);

for i = 1 : no_models
    
    j = i-1;
    
    a11 = hyper_params(i, 1);
    a22 = hyper_params(i, 2);
    
    b11 = hyper_params(i, 3);
    b22 = hyper_params(i, 4);
    
    
    A = makesysd_a(a11, a22, b11, b22, Ts);
    blkA(j*no_states+1 : (j+1)*no_states, j*no_states+1 : (j+1)*no_states) = A.A;
end

blkA_sparse = sparse(blkA);