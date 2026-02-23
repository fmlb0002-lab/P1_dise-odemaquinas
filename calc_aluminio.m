function d = calc_aluminio(T, geom, calib)
% geom: struct con b_mm, e_mm, Lm_mm (solo para epsilon_M)
% calib: struct con K, N0 (puntos para cero)

d.t    = T.("Tiempo");
d.disp = T.("Desplazamiento");
d.F_kN = T.("Fuerza");
d.V_L  = T.("Galga Longitudinal (5)");
d.V_T  = T.("Galga Transversal (6)");

% area inicial
A = geom.b_mm * geom.e_mm;     % mm^2
d.sigma_MPa = 1000*d.F_kN./A;  % MPa (correcto)

% Cero de galgas: fijo si existe, si no calcula
if isfield(calib,'V0_L') && isfield(calib,'V0_T') && ~isempty(calib.V0_L) && ~isempty(calib.V0_T)
    d.V0_L = calib.V0_L;
    d.V0_T = calib.V0_T;
else
    N0 = min(calib.N0, numel(d.V_L));
    d.V0_L = mean(d.V_L(1:N0), "omitnan");
    d.V0_T = mean(d.V_T(1:N0), "omitnan");
end

% deformaciones galgas
d.eps_L = (d.V_L - d.V0_L)*calib.K;
d.eps_T = (d.V_T - d.V0_T)*calib.K;

% deformación "máquina" (solo comparación)
d.eps_M = d.disp./geom.Lm_mm;
end
