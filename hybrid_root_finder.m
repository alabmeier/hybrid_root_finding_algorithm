function [root, info] = hybrid_root_finder(func, Int, params)
% Hybrid Root-Finding Algorithm
%
% Author: Ava Labmeier
%
% Description:
% Combines Inverse Quadratic Interpolation (IQI) and
% Bisection to locate roots of nonlinear functions.
%
% Inputs:
%   func   - function handle
%   Int    - structure containing interval endpoints
%   params - convergence tolerances
%
% Outputs:
%   root   - estimated root
%   info   - convergence information

if ~isfield(params,'max_iter'), params.max_iter = 200; end

a = Int.a; b = Int.b;
fa = func(a); fb = func(b);

if fa*fb > 0
    error('Initial interval must bracket a root: f(a)*f(b) <= 0.');
end
if a == b
    error('Values for a and b must not be equal.');
end
if a > b
    [a,b] = deal(b,a);
    [fa,fb] = deal(fb,fa);
end

% If an endpoint is already a root (within func_tol)
if abs(fa) <= params.func_tol
    root = a;
    info.flag = 0; info.iter = 0; info.fval = fa;
    return
elseif abs(fb) <= params.func_tol
    root = b;
    info.flag = 0; info.iter = 0; info.fval = fb;
    return
end

% Initialize c as midpoint and its value
c = (a + b) / 2;
fc = func(c);

iter = 0;
IQIrepeat = 0;
root = NaN;
info.flag = 1;

while (abs(b - a) > params.root_tol) && (iter < params.max_iter)
    % Try inverse quadratic interpolation (need three distinct x)
    tryIQI = true;
    x0 = a; x1 = b; x2 = c;
    f0 = fa; f1 = fb; f2 = fc;
    % Check distinct f values to avoid division by zero
    if (f0 == f1) || (f0 == f2) || (f1 == f2)
        tryIQI = false;
    end

    if tryIQI
        % compute IQI candidate with safe denominators
        denom0 = (f0 - f1)*(f0 - f2);
        denom1 = (f1 - f0)*(f1 - f2);
        denom2 = (f2 - f0)*(f2 - f1);
        if denom0 == 0 || denom1 == 0 || denom2 == 0
            tryIQI = false;
        else
            x3 = (x0*f1*f2)/denom0 + (x1*f0*f2)/denom1 + (x2*f0*f1)/denom2;
        end
    end

    % If IQI not suitable or candidate outside bracket, use midpoint
    if ~tryIQI || x3 < a || x3 > b || ~isfinite(x3)
        x3 = (a + b) / 2;
        IQIrepeat = 0;
    else
        IQIrepeat = IQIrepeat + 1;
        if IQIrepeat > 5
            % force bisection occasionally to ensure progress
            x3 = (a + b) / 2;
            IQIrepeat = 0;
        end
    end

    fx3 = func(x3);

    % Stopping on function tolerance
    if abs(fx3) <= params.func_tol
        root = x3;
        info.flag = 0;
        info.iter = iter;
        info.fval = fx3;
        return
    end

    % Update bracketing interval: maintain a and b with opposite signs
    if fa * fx3 < 0
        b = x3; fb = fx3;
    else
        a = x3; fa = fx3;
    end

    % Update c as midpoint and its value for next IQI step
    c = (a + b) / 2;
    fc = func(c);

    iter = iter + 1;
end

% Final check
if iter >= params.max_iter
    info.flag = 1;
    info.iter = iter;
    info.fval = NaN;
    root = NaN;
else
    root = (a + b) / 2;
    info.flag = 0;
    info.iter = iter;
    info.fval = func(root);
end
end
