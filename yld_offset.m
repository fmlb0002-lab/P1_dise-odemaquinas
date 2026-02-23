function out = yld_offset(eps, sig, fitE, eps_off, i_start)
eps = eps(:); sig = sig(:);
E = fitE.m; b = fitE.b;

if nargin < 4 || isempty(eps_off); eps_off = 0.002; end
if nargin < 5 || isempty(i_start); i_start = fitE.idx(end); end
i_start = max(2, i_start);

sig_off = E*(eps - eps_off) + b;
d = sig - sig_off;

k = find(d(i_start-1:end-1) > 0 & d(i_start:end) <= 0, 1, 'first');
if ~isempty(k)
    k = k + i_start - 1;
else
    k = find(d(i_start-1:end-1) < 0 & d(i_start:end) >= 0, 1, 'first');
    if isempty(k)
        error('yld_offset: no se encontró cruce offset-curva.');
    end
    k = k + i_start - 1;
end

eps_y = interp1(d(k-1:k), eps(k-1:k), 0);
sig_y = interp1(eps(k-1:k), sig(k-1:k), eps_y);

out.Sy = sig_y;
out.eps = eps_y;
out.sig_off = sig_off;
out.k = k;
end
