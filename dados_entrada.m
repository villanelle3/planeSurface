function [nx, ny, Ly, Re, rho, umax, ucf, CFL, tf, tsave, tsavec, tempoE, ...
         start, nstep_restart, filtro, ci, ab, eb, tipo, BIF, nit, nfi, ...
         n_threads, d] = dados_entrada()

    nx = 1024;		    % nx
    ny = 256;		    % ny
    Ly = 0.38d0;		% Ly dimensional
    Re = 1.0d4;		    % Re
    rho = 1.0d0;		% rho
    umax = 1.0d0;		% umax
    ucf = 0.3d0;		% ucf
    CFL = 0.5d0;		% CFL
    tf = 10.0d0;		% tf - tempo final de simulacao
    tsave = 0.005d0;	% tsave
    tsavec = 1400;		% tsavec campos do escoamento
    tempoE = 200.0d0;	% tempo para comecar as estatisticas
    start = 0;			% start
    nstep_restart = 2;	% nstep_restart
    filtro = 0;			% tipo de filtro
    ci = 5;			    % condicao inicial
    ab = 1.00d0;		% alpha_buffer
    eb = 3;			    % exponecial_buffer
    tipo = 'MDF';		% MDF ou MFV ou BF
    BIF = 0;			% BIF - 1 se variaveis bifasicas; 0 se nao bifasico
    nit = 1000;		    % nit - numero de iteracoes do MDF
    nfi = 0;			% nfi - numero de fronteiras imersas
    n_threads = 1;		% n_threads - numero de threads
    d = 1;			    % diametro

end



