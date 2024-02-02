function [up, vp, pp, wwp, divV] = pos_processamento(fx, fy, nstep, t, dt, u, v, ufis, vfis, iimin, zp, Lx, Ly, nx, ny, dy, jmax, ni, umax, kx, ky, k2)


    % Posição das sondas
    px1 = floor((iimin - zp - 1) - 1); % sz=2 do buffer
    px2 = iimin;
    px3 = iimin + 1;
    px4 = iimin + floor(0.25 * nx / Lx);
    px5 = iimin + floor(0.50 * nx / Lx);
    px6 = iimin + floor(0.75 * nx / Lx);
    
    py1 = floor(0.5 * Ly / dy);
    py2 = jmax;
    py3 = jmax + floor((5.0 * sqrt(ni * 1.0 / umax) / 2.0) / Ly * ny);
    py4 = floor(0.5 * ny);

	ut = u;
    vt = v;
    utfis = ufis;
    vtfis = vfis;
    
	% call tnlinear() !Calculo do tnl para o calculo da pressao
    [tnlxa, tnlya] = tnlinear(ut, vt, utfis, vtfis, nx, ny, kx, ky);
	
    t3 = t;
    p = zeros(nx,ny); cm = zeros(nx,ny); ww = zeros(nx,ny);
    for i = 1:nx
        for j = 1:ny
            p(i, j) = 1i * (1.0 / k2(i, j)) * ((tnlxa(i, j) - fx(i, j)) * kx(i) + (tnlya(i, j) - fy(i, j)) * ky(j));
            cm(i, j) = 1i * kx(i) * u(i, j) + 1i * ky(j) * v(i, j);
            ww(i, j) = 1i * kx(i) * v(i, j) - 1i * ky(j) * u(i, j);
        end
    end

    p = ifft2(p);
    cm = ifft2(cm);
    ww = ifft2(ww);
    
	up = real(ufis);
    vp = real(vfis);
    pp = real(p);
    wwp = real(ww);
    divV = max(max(abs(real(cm))));

    np = 0;
    L2eu = 0.0;
    L2ev = 0.0;

    for i = 1:nx % CC placa plana
        ua(i, jmax) = 0.0;
        va(i, jmax) = 0.0;
    
        u2 = (ua(i, jmax) - up(i, jmax))^2;
        L2eu = L2eu + u2;
    
        v2 = (va(i, jmax) - vp(i, jmax))^2;
        L2ev = L2ev + v2;
    
        np = np + 1;
    end
    L2eu = sqrt(L2eu / np);
    L2ev = sqrt(L2ev / np);

	fprintf('\n'); % Escreve na tela
    fprintf('********************************\n');
    fprintf('tempo: %f\n', t);
    fprintf('nstep: %d\n', nstep);
    fprintf('dt: %f\n', dt);
    fprintf('divV: %f\n', divV);

end 