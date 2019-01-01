function [k_kern, kstar, kstarstar] = kernels(x_data, x_pred)

%% first part of the k(x, x')
[ndp, ~] = size(x_data);

k_kern = zeros(ndp, ndp);

for j = 1 : ndp
   
    diff_mat = x_data - repmat(x_data(j, :), ndp, 1);
    
    k_kern(:, j) = sum((diff_mat.^2), 2);
    
end

%% second part k(x*, x)

[npp, ~] = size(x_pred);

kstar = zeros(npp, ndp);

for j = 1 : ndp
    
   diff_mat = x_pred - repmat(x_data(j, :), npp, 1);
   
   kstar(:, j) = sum((diff_mat.^2), 2);
    
    
end

%% third part k(x*, x*)

kstarstar = zeros(npp, npp);

for j = 1 : npp
    
    
   diff_mat = x_pred - repmat(x_pred(j, :), npp, 1);
   
   kstarstar(:, j) = sum((diff_mat.^2), 2);
    
    
end

end

