function STATS = mcnemar(x, varargin)
%MCNEMAR McNemar's chi-square test for matched pairs (2x2 table).
%
%   mcnemar(x)
%   STATS = mcnemar(x)
%   mcnemar(x, alpha)
%   STATS = mcnemar(x, alpha)
%
%   Description:
%       MCNEMAR performs McNemar's chi-square test on a 2x2 contingency
%       table arising from matched pairs with a dichotomous outcome under
%       two conditions (e.g., pre vs post, drug vs placebo, yes vs no).
%
%       The input x is a 2×2 table:
%
%                     Condition 2
%                     +      -
%           +       x(1,1)  x(1,2)
%   Cond 1
%           -       x(2,1)  x(2,2)
%
%       The test is based on the discordant pairs:
%           b = x(1,2) (positive in condition 1 only)
%           c = x(2,1) (positive in condition 2 only)
%
%       MCNEMAR computes:
%           - The critical chi-square value at significance level alpha
%           - The McNemar chi-square statistic with Yates' continuity
%             correction (1 degree of freedom)
%           - The p-value
%           - An approximate two-sided power
%
%   Inputs:
%       x     - 2×2 numeric matrix of nonnegative integers, real and finite.
%
%       alpha - (Optional) Significance level for the test.
%               Scalar in the interval (0,1). Default: 0.05
%
%   Outputs:
%       STATS - (Optional) structure with fields:
%                  STATS.chisq      : McNemar chi-square with Yates' correction
%                  STATS.df         : degrees of freedom (df = 1)
%                  STATS.crit       : critical chi-square value at alpha
%                  STATS.pvalue     : p-value (two-sided)
%                  STATS.alpha      : significance level
%                  STATS.Zb         : Z_beta (used for power)
%                  STATS.power      : approximate two-sided power
%                  STATS.N          : total number of pairs
%                  STATS.table      : original 2×2 table x
%                  STATS.discordant : [b c] discordant cell counts
%
%       If no output is requested, the function prints a summary of the
%       results to the Command Window.
%
%   Example:
%       % Drug vs placebo example:
%       %                  Drug
%       %             +         -
%       %       +   101       59
%       % Placebo
%       %       -   121       33
%       x = [101 59; 121 33];
%
%       STATS = mcnemar(x);
%
%   GitHub repository:
%       https://github.com/dnafinder/mcnemar
%
%   Citation:
%       Cardillo G. (2025). mcnemar: McNemar's chi-square test for matched
%       pairs in MATLAB. Available at:
%       https://github.com/dnafinder/mcnemar
%
%   License:
%       This function is distributed under the terms specified in the
%       LICENSE file of the mcnemar repository.
%
%   Author:
%       Giuseppe Cardillo
%       giuseppe.cardillo.75@gmail.com
%
%   Created:
%       2007-01-01 (original concept)
%
%   Updated:
%       2025-11-19 (refactored and documented version)
%
%   Version:
%       1.1.0

% -----------------------------
% Input parsing and validation
% -----------------------------
p = inputParser;

% x: 2x2 table of nonnegative integers
addRequired(p, 'x', @(v) validateattributes(v, ...
    {'numeric'}, ...
    {'real', 'finite', 'integer', 'nonnegative', 'nonnan', 'size', [2 2]}));

% alpha: significance level in (0,1)
addOptional(p, 'alpha', 0.05, @(v) ...
    validateattributes(v, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'nonnan', '>', 0, '<', 1}));

parse(p, x, varargin{:});
x     = p.Results.x;
alpha = p.Results.alpha;

clear p;

% -----------------------------
% Extract discordant pairs
% -----------------------------
% Discordant pairs are b and c (off-diagonal elements).
% We can extract them as:
%   b = x(1,2)
%   c = x(2,1)
ob = diag(fliplr(x));  % ob = [b; c]
b  = ob(1);
c  = ob(2);

% Warn if there are no discordant pairs
if b == 0 && c == 0
    warning('mcnemar:NoDiscordantPairs', ...
        'There are no discordant pairs; McNemar''s test cannot be computed.');
end

% -----------------------------
% McNemar chi-square statistic
% -----------------------------
% With Yates' continuity correction and 1 df:
%   chi^2 = (|b - c| - 1)^2 / (b + c)
if (b + c) == 0
    chisq = NaN;
    warning('mcnemar:DegenerateStatistic', ...
        'b + c = 0. McNemar chi-square statistic is undefined.');
else
    chisq = (abs(b - c) - 1)^2 / (b + c);
end

df = 1;
crit = chi2inv(1 - alpha, df);

% p-value (two-sided) from chi-square distribution with 1 df
if isnan(chisq)
    pvalue = NaN;
else
    pvalue = 1 - chi2cdf(chisq, df);
end

% -----------------------------
% Approximate power calculation
% -----------------------------
% Same approximation formula as used in Liddell's function:
Za = abs(-sqrt(2) * erfcinv(alpha));

% Total sample size:
N = sum(x(:));

% Proportion based on minimum discordant cell:
if N == 0
    p = NaN;
else
    p = min(ob ./ N);
end

% Ratio of discordant counts (larger / smaller)
if min(ob) == 0
    pp = Inf;
else
    pp = max(b/c, c/b);
end

if ~isfinite(pp) || isnan(p)
    Zb  = NaN;
    pwr = NaN;
    warning('mcnemar:PowerUndefined', ...
        'Power calculation is undefined when one discordant cell is zero or N = 0.');
else
    num   = abs(sqrt(N * p * (pp - 1)^2) - sqrt(Za^2 * (pp + 1)));
    denom = sqrt(pp + 1 - p * (pp - 1)^2);
    Zb    = num / denom;
    pwr   = (1 - 0.5 * erfc(-Zb / sqrt(2))) * 2;  % two-sided power approximation
end

% -----------------------------
% Display results (if no output)
% -----------------------------
if nargout == 0
    fprintf('Critical value at %0.0f%% significance level = %0.4f\n', ...
        (1 - alpha) * 100, crit);
    fprintf('McNemar chi-square (with Yates'' correction) = %0.6f    p = %0.6f\n', ...
        chisq, pvalue);
    fprintf('alpha = %0.4f  Zb = %0.4f  Power (2-tails) = %0.4f\n', ...
        alpha, Zb, pwr);
end

% -----------------------------
% Build output structure (if requested)
% -----------------------------
if nargout > 0
    STATS.chisq      = chisq;
    STATS.df         = df;
    STATS.crit       = crit;
    STATS.pvalue     = pvalue;
    STATS.alpha      = alpha;
    STATS.Zb         = Zb;
    STATS.power      = pwr;
    STATS.N          = N;
    STATS.table      = x;
    STATS.discordant = [b, c];
end

end
