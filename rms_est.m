function [index, rms_vals] = rms_est(y_mat, y_data)

diff_mat = y_mat-y_data;

shape_diff_mat = size(diff_mat);

diff_mat = diff_mat.^2;
sq_mat = sum(diff_mat, 2)/shape_diff_mat(2);

[~, index] = min(sq_mat, [], 1);
rms_vals = sqrt(sq_mat);
