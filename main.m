clear; clc;
% -------------------------------------------------------------------------
[nx, ny, Ly, Re, rho, umax, ucf, CFL, tf, tsave, tsavec, tempoE, start,...
nstep_restart, filtro, t_ci, ab, eb, tipo, BIF, nit, nfi, n_threads, d] = dados_entrada();

if BIF == 1, warning("Bifásico ainda não implementado!"); end
% -------------------------------------------------------------------------

% Viscosidade imposta
ni = (umax)/Re; 
% -------------------------------------------------------------------------

% Parâmetros globais
[Lx, dx, dy, x, y, kx, ky, k2, Sx, Sy, gx, gy, tsup, Eo, M, lambda_r, ...
lambda_m, La, teta, fx, fy, rr] = parametros_globais(Ly, nx, ny, d);
% -------------------------------------------------------------------------

% Geracao da geometria a ser imersa
[imin, jmax, iimin, iimax, jjmax, zp] = gera_geometria(nx, ny, Ly);
% -------------------------------------------------------------------------

% Condicao inicial
[u, v, p] = ci(t_ci, nx, ny, Sx, Sy, ni, 0, umax, ucf, rr, teta, jmax, jjmax);


