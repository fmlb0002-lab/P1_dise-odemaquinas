function fit = ajuste_lin(x, y, i0, i1)
%AJUSTE_LIN Ajuste lineal y = m*x + b en el rango i0:i1 (índices)
% Devuelve struct: fit.m, fit.b, fit.R2, fit.x, fit.y, fit.yfit, fit.idx

x = x(:); 
y = y(:);
N = numel(x);

if numel(y) ~= N
    error('ajuste_lin: x e y deben tener la misma longitud.');
end
if i0 < 1 || i1 > N || i0 >= i1
    error('ajuste_lin: rango inválido (i0=%d, i1=%d, N=%d).', i0, i1, N);
end

idx = (i0:i1).';
xx  = x(idx);
yy  = y(idx);

% Limpieza mínima
ok = isfinite(xx) & isfinite(yy);
xx = xx(ok);
yy = yy(ok);
idx = idx(ok);

p    = polyfit(xx, yy, 1);
m    = p(1);
b    = p(2);
yfit = polyval(p, xx);

SSres = sum((yy - yfit).^2);
SStot = sum((yy - mean(yy)).^2);
R2    = 1 - SSres/SStot;

fit.m = m;
fit.b = b;
fit.R2 = R2;
fit.x = xx;
fit.y = yy;
fit.yfit = yfit;
fit.idx = idx;
end
