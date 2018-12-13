function [index, rms_vals] = rms_est(y_mat, y_data)

[y_mat_down, ~] = size(y_mat);
[y_down, ~] = size(y_data);

no_reps = y_mat_down / y_down;
y_reshape = repmat(y_data, no_reps, 1);

diff_mat = y_mat - y_reshape;

shape_diff_mat = size(diff_mat);

diff_mat = diff_mat.^2;
sq_mat = sum(diff_mat, 2)/shape_diff_mat(2);

[~, index] = min(sq_mat, [], 1);
rms_vals = sqrt(sq_mat);
