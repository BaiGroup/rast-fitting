function err = relative_error_safe(x_predicted, x_true)
err = x_predicted - x_true;
idx_imperfect = err ~= 0;
idx_zero_true = (x_true == 0) & idx_imperfect;
idx_nonzero_true = (x_true ~= 0) & idx_imperfect;
err(idx_zero_true) = err(idx_zero_true) ./ x_predicted(idx_zero_true);
err(idx_nonzero_true) = err(idx_nonzero_true) ./ x_true(idx_nonzero_true);
end
