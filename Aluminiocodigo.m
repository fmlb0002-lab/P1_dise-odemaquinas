clearvars; close all; clc

% ===== Config aluminio =====
fname = "aluminio_1.csv";  %<<-- NOMBRE DEL ARCHIVO DE DATOS

geom.b_mm = 12.6;   %(mm)
geom.e_mm = 2.1;    %(mm)
geom.Lm_mm = 75;    %(mm)  SOLO para eps_M (comparación)

calib.K  = 25600e-6;  % (mm/mm)/V  (valor anotado en laboratorio)
calib.V0_L = 1.2431;   % fijo
calib.V0_T = 1.2232;   % fijo

%====================================
%TRAMO ELASTICO, HAY QUE CONFIGURARLO VIENDO COMO SALE R^2
i0 = 5;
i1 = 35;
% ===================================

T = cargar_csv_ensayo(fname);
d = calc_aluminio(T, geom, calib);

% ===== Cálculo del módulo elástico (E) por ajuste lineal en el tramo [i0, i1] =====
fitE = ajuste_lin(d.eps_L, d.sigma_MPa, i0, i1);

E_MPa = fitE.m;          % porque sigma está en MPa y epsilon es adimensional
E_GPa = E_MPa/1000;

fprintf("E = %.3f GPa\n", E_GPa);
fprintf("Ajuste: sigma = %.6g*eps_L + %.6g  (R^2 = %.5f)\n", fitE.m, fitE.b, fitE.R2);
fprintf("Tramo usado: i=[%d,%d], eps_L=[%.3e, %.3e]\n", i0, i1, d.eps_L(i0), d.eps_L(i1));

