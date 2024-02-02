function [up, vp, pp, wwp, divV] = pos_processamento(fx, fy, nstep, t, dt, u, v, ufis, vfis, iimin, zp, Lx, Ly, nx, ny, dy, jmax, ni, umax, kx, ky, k2)
    % Variable initialization
    p = zeros(nx, ny);
    cm = zeros(nx, ny);
    ww = zeros(nx, ny);
    L2eu = 0;
    L2ev = 0;

    % Position of probes
    px1 = int32((iimin - zp - 1) - 1); % sz=2 do buffer
    px2 = iimin;
    px3 = iimin + 1;
    px4 = iimin + int32(0.25 / Lx * nx);
    px5 = iimin + int32(0.50 / Lx * nx);
    px6 = iimin + int32(0.75 / Lx * nx);

    py1 = int32(0.5 * Ly / dy);
    py2 = jmax;
    py3 = jmax + int32((5.0 * sqrt(ni / umax) / 2.0) / Ly * ny);
    py4 = int32(0.5 * ny);

    % FFT
    ut = u;
    vt = v;
    utfis = ufis;
    vtfis = vfis;
    [tnlxa, tnlya] = tnlinear(ut, vt, utfis, vtfis, nx, ny, kx, ky);

    % Pressure calculation
    for i = 1:nx
        for j = 1:ny
            p(i, j) = 1i * 1.0 / k2(i, j) * ((tnlxa(i, j) - fx(i, j)) * kx(i) + (tnlya(i, j) - fy(i, j)) * ky(j)); % -fontex(i, j) -fontey(i, j)
            cm(i, j) = 1i * kx(i) * u(i, j) + 1i * ky(j) * v(i, j);
            ww(i, j) = 1i * kx(i) * v(i, j) - 1i * ky(j) * u(i, j);
        end
    end

    % Transform back
    p = ifft2(p);
    cm = ifft2(cm);
    ww = ifft2(ww);
    up = real(ufis);
    vp = real(vfis);
    pp = real(p);
    wwp = real(ww);
    divV = max(max(abs(real(cm))));

    % L2 Norm calculations
    np = 0;
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

    % blasius(); % Call to Blasius subroutine

    % Output
    disp(' '); % Writes to the console
    disp('********************************');
    disp(['tempo:', num2str(t)]);
    disp(['nstep:', num2str(nstep)]);
    disp(['dt:', num2str(dt)]);
    disp(['divV:', num2str(divV)]);
end
