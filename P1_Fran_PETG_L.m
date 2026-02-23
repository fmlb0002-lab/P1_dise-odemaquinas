clearvars; close all; clc

% ===== Config petG =====
fname1 = "PETG CF Trans_1.csv";  %<<-- NOMBRE DEL ARCHIVO DE DATOS
fname2 = "PETGCF_long_1.csv"; 

geom.b_mm = 13.33;   %(mm)
geom.e_mm = 4;    %(mm)
geom.Lm_mm = 57;    %(mm)  SOLO para eps_M (comparación)

calib.K  = 25600e-6;  % (mm/mm)/V  (valor anotado en laboratorio)

calib.V0_L1 = 1.2345;   % fijo PROB 1
calib.V0_T1 = 1.2299;   % fijo

calib.V0_L2 = 1.2351;   % fijo PROB2
calib.V0_T2 = 1.2298;   % fijo
%====================================
%TRAMO ELASTICO, HAY QUE CONFIGURARLO VIENDO COMO SALE R^2
i1_0 = 10;
i1_f = 60;
i2_0 = 10;
i2_f = 60;
% ===================================


%Cargo los datos de ambas provetas 


% --- Probeta 1 (Trans) ---
calib1 = calib;
calib1.V0_L = calib.V0_L1;
calib1.V0_T = calib.V0_T1;
T1 = cargar_csv_ensayo(fname1);
d1 = calc_petg(T1, geom, calib1);

% --- Probeta 1 (Long) ---

calib2 = calib;
calib2.V0_L = calib.V0_L2;
calib2.V0_T = calib.V0_T2;
T2 = cargar_csv_ensayo(fname2);
d2 = calc_petg(T2, geom, calib2);

%OJO!!! en los nombres a la hora de conectar el longitudinal y transversal
%se conectaron de forma incorrecta, por lo que vamos a hacer el cambio
aux1 = d1.eps_L;
d1.eps_L= d1.eps_T; 
d1.eps_T = aux1;  

aux2 = d2.eps_L;
d2.eps_L = d2.eps_T; 
d2.eps_T = aux2;

%Hacemos los ajustes lineales de ambas probetas

fitE1 = ajuste_lin(d1.eps_L,d1.sigma_MPa,i1_0,i1_f);
fitE2 = ajuste_lin(d2.eps_L, d2.sigma_MPa, i2_0, i2_f);



%Formamos las rectas para la representacion posterior
x1 = d1.eps_L(i1_0:i1_f);  y1 = d1.sigma_MPa(i1_0:i1_f);
x2 = d2.eps_L(i2_0:i2_f);  y2 = d2.sigma_MPa(i2_0:i2_f);
y1_fit = fitE1.m.*x1 + fitE1.b;
y2_fit = fitE2.m.*x2 + fitE2.b;

%Calculo de las magnitudes mas interesantes.
% ===== 1) MODULO DE YOUNG (pendiente del ajuste) =====
E1_MPa = fitE1.m;          % MPa
E2_MPa = fitE2.m;          % MPa
E1_GPa = E1_MPa/1000;      % GPa
E2_GPa = E2_MPa/1000;      % GPa

fprintf('\n=== MODULO DE YOUNG ===\n');
fprintf('E1 (Trans) = %.2f MPa  = %.3f GPa   (R2=%.5f)\n', E1_MPa, E1_GPa, fitE1.R2);
fprintf('E2 (Long)  = %.2f MPa  = %.3f GPa   (R2=%.5f)\n', E2_MPa, E2_GPa, fitE2.R2);



% =====2)  Sy por OFFSET 0.2% (sigma_0.2) =====
eps_off = 0.002;

% Toe compensation (cero corregido) + recta elástica pasando por (0,0)
eps0_1 = -fitE1.b/fitE1.m;
eps0_2 = -fitE2.b/fitE2.m;

eps1 = d1.eps_L - eps0_1;  sig1 = d1.sigma_MPa;
eps2 = d2.eps_L - eps0_2;  sig2 = d2.sigma_MPa;

fitE1c = fitE1; fitE1c.b = 0;
fitE2c = fitE2; fitE2c.b = 0;

% Offset yield
yld1 = yld_offset(eps1, sig1, fitE1c, eps_off, i1_f);
yld2 = yld_offset(eps2, sig2, fitE2c, eps_off, i2_f);

fprintf('\n=== 0.2%% OFFSET YIELD ===\n');
fprintf('Transversal:  sigma_0.2 = %.2f MPa (eps=%.6f)\n', yld1.Sy, yld1.eps);
fprintf('Longitudinal: sigma_0.2 = %.2f MPa (eps=%.6f)\n', yld2.Sy, yld2.eps);

% =====3)  ´Tensión máxima Su =====
Su1 = max(d1.sigma_MPa, [], 'omitnan');
Su2 = max(d2.sigma_MPa, [], 'omitnan');

fprintf("Su Trans = %.2f MPa\n", Su1);
fprintf("Su Long  = %.2f MPa\n", Su2);

% =====4)  Tension de rotura Sf =====
% ===== Sf: última tensión "válida" con sigma > 10 MPa =====
thr = 10; % MPa discrimino a partir de este valor para excluir valores que no sean de rotua, en el caso de que coincida se debera cambiar este valor
iSf1 = find(d1.sigma_MPa > thr, 1, 'last');
iSf2 = find(d2.sigma_MPa > thr, 1, 'last');
Sf1 = d1.sigma_MPa(iSf1);
Sf2 = d2.sigma_MPa(iSf2);
epsf1 = d1.eps_L(iSf1);
epsf2 = d2.eps_L(iSf2);

fprintf("Sf Trans = %.2f MPa (eps=%.6f)\n", Sf1, epsf1);
fprintf("Sf Long  = %.2f MPa (eps=%.6f)\n", Sf2, epsf2);

% ===== 5) Poisson ν con ajuste_lin en tramo elástico =====

% --- Probeta 1 ---
x1 = d1.eps_L(i1_0:i1_f);
y1 = d1.eps_T(i1_0:i1_f);

x1 = x1 - x1(1);
y1 = y1 - y1(1);

fitNu1 = ajuste_lin(x1, y1, 1, numel(x1));
nu1 = -fitNu1.m;

% --- Probeta 2 ---
x2 = d2.eps_L(i2_0:i2_f);
y2 = d2.eps_T(i2_0:i2_f);

x2 = x2 - x2(1);
y2 = y2 - y2(1);

fitNu2 = ajuste_lin(x2, y2, 1, numel(x2));
nu2 = -fitNu2.m;

fprintf('\n=== POISSON ===\n');
fprintf('nu1 (Trans) = %.4f   (R2=%.5f)\n', nu1, fitNu1.R2);
fprintf('nu2 (Long)  = %.4f   (R2=%.5f)\n', nu2, fitNu2.R2);


% ===== 6) Integración numérica hasta eps_f =====
epsf_T = 0.013010;   % Trans 
epsf_L = 0.012946;   % Long  

% Transversal
mask1 = d1.eps_L <= epsf_T;
U1 = trapz(d1.eps_L(mask1), d1.sigma_MPa(mask1));   % MPa = MJ/m^3

% Longitudinal
mask2 = d2.eps_L <= epsf_L;
U2 = trapz(d2.eps_L(mask2), d2.sigma_MPa(mask2));

fprintf("U Trans (hasta epsf)= %.2f MJ/m^3\n", U1);
fprintf("U Long  (hasta epsf)= %.2f MJ/m^3\n", U2);


figure;
plot(d1.eps_M, d1.sigma_MPa, 'LineWidth', 1.1); hold on;
plot(d2.eps_M, d2.sigma_MPa, 'LineWidth', 1.1); hold off;
grid on;
xlabel('\epsilon_M = \Delta L / L_m');
ylabel('\sigma (MPa)');
title('PETG-CF: Transversal vs Longitudinal');
legend('Transversal','Longitudinal','Location','best');


%Representamos DATOS + AJUSTE en la misma gráfica

figure;
plot(d1.eps_L, d1.sigma_MPa, '.', 'MarkerSize', 10); hold on;
plot(d2.eps_L, d2.sigma_MPa, '.', 'MarkerSize', 10);
plot(x1, y1_fit, '-', 'LineWidth', 4);
plot(x2, y2_fit, '-', 'LineWidth', 4);
grid on;
xlabel('\epsilon_L (galga longitudinal)');
ylabel('\sigma (MPa)');
title('PETG-CF: Transversal vs Longitudinal (puntos + ajuste lineal)');
legend( ...
  'Transversal (datos)', ...
  'Longitudinal (datos)', ...
  sprintf('Transversal: \\sigma = %.4g\\epsilon + %.4g (R^2=%.5f)', fitE1.m, fitE1.b, fitE1.R2), ...
  sprintf('Longitudinal: \\sigma = %.4g\\epsilon + %.4g (R^2=%.5f)', fitE2.m, fitE2.b, fitE2.R2), ...
  'Location','best');
hold off

%============================================================================
%/////////////////       REPRESENTACION DE FIGURAS       ///////////////// 
%============================================================================


clearvars; close all; clc

%% =========================
%  PETG-CF (ASTM D638) — Transversal vs Longitudinal
%  Script ordenado: 1) carga/cálculo  2) magnitudes  3) figuras
%% =========================

%% 0) Configuración
fnameT = "PETG CF Trans_1.csv";
fnameL = "PETGCF_long_1.csv";

geom.b_mm  = 13.33;  % mm
geom.e_mm  = 4.00;   % mm
geom.Lm_mm = 57;     % mm (solo para eps_M)

calib.K = 25600e-6;  % (mm/mm)/V

% Ceros por probeta (fijos)
calibT.V0_L = 1.2345;  calibT.V0_T = 1.2299;  % Trans
calibL.V0_L = 1.2351;  calibL.V0_T = 1.2298;  % Long

% Tramo elástico (índices) para E y ν
iT0 = 10; iTf = 60;
iL0 = 10; iLf = 60;

% Parámetros comunes
eps_off = 0.002;   % 0.2%
thr_sig = 10;      % MPa (umbral para definir "fin" de señal útil)

%% 1) Carga de datos y cálculo de magnitudes base
TT = cargar_csv_ensayo(fnameT);
TL = cargar_csv_ensayo(fnameL);

dT = calc_petg(TT, geom, struct('K',calib.K,'V0_L',calibT.V0_L,'V0_T',calibT.V0_T));
dL = calc_petg(TL, geom, struct('K',calib.K,'V0_L',calibL.V0_L,'V0_T',calibL.V0_T));

% Corrección de cableado: longitudinal/transversal conectados al revés
dT = swap_epsLT(dT);
dL = swap_epsLT(dL);

%% 2) Cálculo de parámetros (sin figuras)

% 2.1) Módulo de Young (E): ajuste lineal en tramo [i0,if]
fitET = ajuste_lin(dT.eps_L, dT.sigma_MPa, iT0, iTf);
fitEL = ajuste_lin(dL.eps_L, dL.sigma_MPa, iL0, iLf);

ET_GPa = fitET.m/1000;
EL_GPa = fitEL.m/1000;

fprintf('\n=== MODULO DE YOUNG ===\n');
fprintf('E (Trans) = %.2f MPa = %.3f GPa (R2=%.5f)\n', fitET.m, ET_GPa, fitET.R2);
fprintf('E (Long)  = %.2f MPa = %.3f GPa (R2=%.5f)\n', fitEL.m, EL_GPa, fitEL.R2);

% 2.2) 0.2% offset yield (toe compensation + recta elástica pasando por (0,0))
%      Toe compensation: desplazar eps para que el ajuste elástico cruce por (0,0)
eps0_T = -fitET.b/fitET.m;
eps0_L = -fitEL.b/fitEL.m;

epsT = dT.eps_L - eps0_T;  sigT = dT.sigma_MPa;
epsL = dL.eps_L - eps0_L;  sigL = dL.sigma_MPa;

fitET0 = fitET; fitET0.b = 0;
fitEL0 = fitEL; fitEL0.b = 0;

i_startT = find(epsT >= eps_off, 1, 'first');
i_startL = find(epsL >= eps_off, 1, 'first');

yldT = yld_offset(epsT, sigT, fitET0, eps_off, i_startT);
yldL = yld_offset(epsL, sigL, fitEL0, eps_off, i_startL);

fprintf('\n=== 0.2%% OFFSET YIELD ===\n');
fprintf('Transversal:  sigma_0.2 = %.2f MPa (eps=%.6f)\n', yldT.Sy, yldT.eps);
fprintf('Longitudinal: sigma_0.2 = %.2f MPa (eps=%.6f)\n', yldL.Sy, yldL.eps);

% 2.3) Tensión máxima (ISO: Rm, aquí lo llamas Su)
[SuT, iRmT] = max(dT.sigma_MPa, [], 'omitnan');
[SuL, iRmL] = max(dL.sigma_MPa, [], 'omitnan');

fprintf('\n=== MAXIMO ===\n');
fprintf('Su Trans = %.2f MPa\n', SuT);
fprintf('Su Long  = %.2f MPa\n', SuL);

% 2.4) Punto final (rotura): último punto tras el máximo con sigma > thr_sig
iEndT_rel = find(dT.sigma_MPa(iRmT:end) > thr_sig, 1, 'last');
iEndL_rel = find(dL.sigma_MPa(iRmL:end) > thr_sig, 1, 'last');

iSfT = iRmT + iEndT_rel - 1;
iSfL = iRmL + iEndL_rel - 1;

SfT   = dT.sigma_MPa(iSfT);   epsfT = dT.eps_L(iSfT);
SfL   = dL.sigma_MPa(iSfL);   epsfL = dL.eps_L(iSfL);

fprintf('\n=== PUNTO FINAL (rotura operativa) ===\n');
fprintf('Sf Trans = %.2f MPa (eps=%.6f)\n', SfT, epsfT);
fprintf('Sf Long  = %.2f MPa (eps=%.6f)\n', SfL, epsfL);

% 2.5) Poisson: ajuste de -eps_T vs eps_L en el mismo tramo elástico
fitNuT = ajuste_lin(dT.eps_L, -dT.eps_T, iT0, iTf);
fitNuL = ajuste_lin(dL.eps_L, -dL.eps_T, iL0, iLf);

nuT = fitNuT.m;
nuL = fitNuL.m;

fprintf('\n=== POISSON ===\n');
fprintf('nu (Trans) = %.4f (R2=%.5f)\n', nuT, fitNuT.R2);
fprintf('nu (Long)  = %.4f (R2=%.5f)\n', nuL, fitNuL.R2);

% 2.6) Energía específica U = ∫σ dε hasta epsf (con eps_L)
maskT = dT.eps_L <= epsfT;
maskL = dL.eps_L <= epsfL;

UT = trapz(dT.eps_L(maskT), dT.sigma_MPa(maskT));   % MPa = MJ/m^3
UL = trapz(dL.eps_L(maskL), dL.sigma_MPa(maskL));

fprintf('\n=== ENERGIA ===\n');
fprintf('U Trans (hasta epsf)= %.2f MJ/m^3\n', UT);
fprintf('U Long  (hasta epsf)= %.2f MJ/m^3\n', UL);

%% 3) Figuras (aquí ya no se recalcula nada)

% FIG P1: σ-ε superpuestas con puntos característicos
figure('Color','w')
plot(dT.eps_L, dT.sigma_MPa, '.', 'MarkerSize', 9); hold on
plot(dL.eps_L, dL.sigma_MPa, '.', 'MarkerSize', 9);

plot(dT.eps_L(iRmT), SuT, 's', 'MarkerSize', 8, 'LineWidth', 1.2)
plot(dL.eps_L(iRmL), SuL, 's', 'MarkerSize', 8, 'LineWidth', 1.2)

plot(yldT.eps + eps0_T, yldT.Sy, 'o', 'MarkerSize', 8, 'LineWidth', 1.2) % reubico a eps original
plot(yldL.eps + eps0_L, yldL.Sy, 'o', 'MarkerSize', 8, 'LineWidth', 1.2)

plot(epsfT, SfT, '^', 'MarkerSize', 8, 'LineWidth', 1.2)
plot(epsfL, SfL, '^', 'MarkerSize', 8, 'LineWidth', 1.2)

grid on
xlabel('\epsilon')
ylabel('\sigma (MPa)')
title('PETG-CF: curvas \sigma-\epsilon (Transversal vs Longitudinal)')
legend('Datos Trans','Datos Long', ...
       'Su Trans','Su Long', ...
       '\sigma_{0.2} Trans','\sigma_{0.2} Long', ...
       '(\epsilon_f,\sigma_f) Trans','(\epsilon_f,\sigma_f) Long', ...
       'Location','best')
hold off
% exportgraphics(gcf,'FIGURAS/petg_curvas_superpuestas.png','Resolution',300);

% FIG P2: Zoom elástico + ajuste E
figure('Color','w')
subplot(1,2,1)
plot(dT.eps_L, dT.sigma_MPa, '.', 'MarkerSize', 9); hold on
plot(fitET.x, fitET.yfit, '-', 'LineWidth', 1.8)
grid on; xlabel('\epsilon_L'); ylabel('\sigma (MPa)')
title(sprintf('Trans: E=%.3f GPa (R^2=%.5f)', ET_GPa, fitET.R2))
xlim([0 dT.eps_L(iTf)*1.2]); ylim([0 max(dT.sigma_MPa(iT0:iTf))*1.3]); hold off

subplot(1,2,2)
plot(dL.eps_L, dL.sigma_MPa, '.', 'MarkerSize', 9); hold on
plot(fitEL.x, fitEL.yfit, '-', 'LineWidth', 1.8)
grid on; xlabel('\epsilon_L'); ylabel('\sigma (MPa)')
title(sprintf('Long: E=%.3f GPa (R^2=%.5f)', EL_GPa, fitEL.R2))
xlim([0 dL.eps_L(iLf)*1.2]); ylim([0 max(dL.sigma_MPa(iL0:iLf))*1.3]); hold off
% exportgraphics(gcf,'FIGURAS/petg_E_zoom.png','Resolution',300);

% FIG P3: Offset 0.2%
eps_zoom = 0.04;

figure('Color','w')
subplot(1,2,1)
mask = dT.eps_L>=0 & dT.eps_L<=eps_zoom;
plot(dT.eps_L(mask), dT.sigma_MPa(mask), '.', 'MarkerSize', 9); hold on
plot(fitET.x, fitET.yfit, '-', 'LineWidth', 1.4)
% recta offset evaluada en eps corregida -> la pinto sobre eps original sumando eps0_T
plot((dT.eps_L(mask)-eps0_T), yldT.sig_off(mask), '--', 'LineWidth', 1.4)
plot(yldT.eps + eps0_T, yldT.Sy, 'o', 'MarkerSize', 8, 'LineWidth', 1.2)
grid on; xlabel('\epsilon_L'); ylabel('\sigma (MPa)')
title(sprintf('Trans: \\sigma_{0.2}=%.2f MPa', yldT.Sy))
xlim([0 eps_zoom]); ylim([0 max(dT.sigma_MPa(mask))*1.2])
legend('Datos','Recta elástica','Offset 0.2%','\sigma_{0.2}','Location','best')
hold off

subplot(1,2,2)
mask = dL.eps_L>=0 & dL.eps_L<=eps_zoom;
plot(dL.eps_L(mask), dL.sigma_MPa(mask), '.', 'MarkerSize', 9); hold on
plot(fitEL.x, fitEL.yfit, '-', 'LineWidth', 1.4)
plot((dL.eps_L(mask)-eps0_L), yldL.sig_off(mask), '--', 'LineWidth', 1.4)
plot(yldL.eps + eps0_L, yldL.Sy, 'o', 'MarkerSize', 8, 'LineWidth', 1.2)
grid on; xlabel('\epsilon_L'); ylabel('\sigma (MPa)')
title(sprintf('Long: \\sigma_{0.2}=%.2f MPa', yldL.Sy))
xlim([0 eps_zoom]); ylim([0 max(dL.sigma_MPa(mask))*1.2])
legend('Datos','Recta elástica','Offset 0.2%','\sigma_{0.2}','Location','best')
hold off
% exportgraphics(gcf,'FIGURAS/petg_yield_offset.png','Resolution',300);

% FIG P4: Poisson (-eps_T vs eps_L) tramo elástico
figure('Color','w')
subplot(1,2,1)
idx = iT0:iTf;
plot(dT.eps_L(idx), -dT.eps_T(idx), '.', 'MarkerSize', 10); hold on
plot(fitNuT.x, fitNuT.yfit, '-', 'LineWidth', 1.8)
grid on; xlabel('\epsilon_L'); ylabel('-\epsilon_T')
title(sprintf('Trans: \\nu=%.4f (R^2=%.5f)', nuT, fitNuT.R2))
hold off

subplot(1,2,2)
idx = iL0:iLf;
plot(dL.eps_L(idx), -dL.eps_T(idx), '.', 'MarkerSize', 10); hold on
plot(fitNuL.x, fitNuL.yfit, '-', 'LineWidth', 1.8)
grid on; xlabel('\epsilon_L'); ylabel('-\epsilon_T')
title(sprintf('Long: \\nu=%.4f (R^2=%.5f)', nuL, fitNuL.R2))
hold off
% exportgraphics(gcf,'FIGURAS/petg_nu.png','Resolution',300);

% FIG P5: Barras resumen (E, σ0.2, Su, epsf, U)
valsT = [ET_GPa, yldT.Sy, SuT];
valsL = [EL_GPa, yldL.Sy, SuL];
labels = {'E (GPa)','\sigma_{0.2} (MPa)','Su (MPa)'};

figure('Color','w')
bar([valsT; valsL]')
grid on
set(gca,'XTickLabel',labels)
xtickangle(25)
ylabel('Valor')
title('PETG-CF: comparación de propiedades (Trans vs Long)')
legend('Transversal','Longitudinal','Location','best')
% exportgraphics(gcf,'FIGURAS/petg_barras_resumen.png','Resolution',300);

%% --------- función local: intercambio eps_L / eps_T ----------
function d = swap_epsLT(d)
    aux = d.eps_L;
    d.eps_L = d.eps_T;
    d.eps_T = aux;
end
