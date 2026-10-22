%fits a line to the error data in log space to get p (order) and k (constant)
%from the scaling relation error = k*(x)^p, same fit we used in Mod 1
%(x is the step size h or the number of rate function calls)
%INPUTS:
%   x_regression: step sizes (or number of rate function calls)
%   y_regression: truncation error for each one
%OUTPUTS:
%   p: order of the scaling (slope of the fit)
%   k: scaling constant (intercept of the fit, undone from log)
function [p,k] = generate_error_fit(x_regression,y_regression)
    x_regression = x_regression(:);
    y_regression = y_regression(:);
    valid = isfinite(x_regression) & isfinite(y_regression) & ...
            x_regression > 0 & y_regression > 0;
    x_regression = x_regression(valid);
    y_regression = y_regression(valid);

    if length(x_regression) < 3 || length(unique(x_regression)) < 2
        error('Not enough error pairs to fit p and k');
    end

    Y = log(y_regression);
    X1 = log(x_regression);
    X2 = ones(length(X1),1);

    coeff_vec = [X1,X2]\Y;

    p = coeff_vec(1);
    k = exp(coeff_vec(2));
end
