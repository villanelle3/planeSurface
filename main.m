clear; clc;
% -------------------------------------------------------------------------
% Load input data
[nx, ny, Ly, Re, rho, umax, ucf, cfl, tf, tsave, nsavec, tempoE, start, ...
    nstep_restart, filtro, t_ci, alpha_buffer, eb, tipo, BIF, nit, nfi, ...
    n_threads, d] = dados_entrada();

if BIF == 1
    warning("Bifásico ainda não implementado!"); 
end
% -------------------------------------------------------------------------

% Imposed viscosity
ni = umax / Re; 
% -------------------------------------------------------------------------

% Global parameters
[Lx, dx, dy, x, y, kx, ky, k2, Sx, Sy, gx, gy, tsup, Eo, M, lambda_r, ...
    lambda_m, La, teta, fx, fy, rr] = parametros_globais(Ly, nx, ny, d);
% -------------------------------------------------------------------------

% Generate immersed geometry
[imin, jmax, iimin, iimax, jjmax, zp] = gera_geometria(nx, ny, Ly);
% -------------------------------------------------------------------------

% Initial conditions
[u, v, p] = ci(t_ci, nx, ny, Sx, Sy, ni, 0, umax, ucf, rr, teta, jmax, jjmax);
% -------------------------------------------------------------------------

% Calculate buffer zone
[phi, Qtx, Qty, Qtxb, Qtyb, b1] = buffer(nx, ny, kx, dx, imin, iimin, ...
    alpha_buffer, x, y, jmax, jjmax, umax, ni);
% -------------------------------------------------------------------------

ufis = real(u);
vfis = real(v);
umax = max(max(abs(ufis)));
vmax = max(max(abs(vfis)));

dtd = 2.0 / ni * (1.0 / (dx^2) + 1.0 / (dy^2))^-1;    % dt difusivo
dta = min(dx / (umax + eps / 2), dy / (vmax + eps / 2));    % dt advectivo
dt = cfl * min(dta, dtd); % Discretizacao do tempo

tmin = 0.0;
t1 = 0.0;
t2 = 0.0;
t = 0.0;
nstep = 0;
nstepc = 1;
% -------------------------------------------------------------------------

u = fft2(u);
v = fft2(v);
um = 0.0;
vm = 0.0;
r11 = 0.0;
r12 = 0.0;
r22 = 0.0;
nn = 0;
nn1 = 0;
nn2 = 0;
ie = 0; % int(te/dt) !iteracao que comeca as estatisticas
% -------------------------------------------------------------------------

fprintf('Re=%f\n', Re);
fprintf('ni=%f\n', ni);
fprintf('dx=%f\n', dx);
fprintf('dy=%f\n', dy);
% fprintf('ds=%f\n', front(1).ds);
% fprintf('nl=%f\n', front(1).nl);
% fprintf('MA=%f\n', front(1).ma);
fprintf('dt=%f\n', dt);
% -------------------------------------------------------------------------

while t <= tf + dt % Início do avanço temporal
    dta = min(dx / max(abs(real(ufis(:)))), dy / max(abs(real(vfis(:))))); % dt advectivo
    dt = cfl * min(dta, dtd); % Discretização do tempo

    [up, vp, pp, wwp, divV] = pos_processamento(fx, fy, nstep, t, dt, u, ...
        v, ufis, vfis, iimin, zp, Lx, Ly, nx, ny, dy, jmax, ni, umax, kx, ...
        ky, k2); % Subrotina que gera os arquivos de saída

    % Esquema de avanço temporal (subrotina avanco_temporal.f90)
    [ufis, vfis, u, v] = RK46DF(dt, ni, nx, ny, ufis, vfis, kx, ky, k2, u, v);
    [fx, fy, u, v, ufis, vfis] = forcaDF(u, v, ufis, vfis, nx,...
        ny, phi, Qtx, Qty, kx, ky, k2, nit, jmax, jjmax, dt, zp, rho,...
        nstep, nsavec, iimin, iimax, dy, umax, ni);
    [u, v] = projecao(u, v, nx, ny, kx, ky, k2);

    nstep = nstep + 1; % Incremento de nstep
    t = t + dt; % nstep * dt; % Incremento do tempo de dt
end % Fim do avanço temporal
