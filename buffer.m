function [phi, Qtx, Qty, Qtxb, Qtyb, b1] = buffer(nx, ny, kx, dx, imin, iimin, alpha_buffer, x, y, jmax, jjmax, umax, ni)
    sz = imin;
    drise = 0.6 * (x(iimin) - x(sz));
    dfall = 0.2 * (x(iimin) - x(sz));
    phi = zeros(nx, ny);

    for i = sz:iimin
        for j = jmax + 1:jjmax - 1
            S1 = funcS((x(i) - x(sz)) / drise);
            S2 = funcS((x(i) - x(iimin)) / dfall + 1);
            phi(i, j) = alpha_buffer * (S1 - S2);
        end
    end

    Qtx = zeros(nx, ny);
    Qty = zeros(nx, ny);

    for i = 1:iimin
        xf = x(i) - x(iimin) + x(nx) + 0.5 + dx;
        xi = 0.5 - ((iimin) - i) * dx;
        S1 = funcS((x(i) - x(sz)) / ((x(iimin) - x(sz))));
        xr(i) = xf + (xi - xf) * S1;
    end

    for i = iimin + 1:nx
        xr(i) = x(i) - x(iimin) + 0.5;
    end

    dxrdx = xr;
    for i = 1:nx
        dxrdx(i) = 1i * kx(i) * dxrdx(i);
    end
    b1 = real(ifft2(dxrdx));

    for i = 1:nx
        [uLoc, vLoc] = perfilLocal(xr(i), ny, y, umax, ni);
        for j = jmax:jjmax
            Qtx(i, j) = uLoc(j);
            Qty(i, j) = dxrdx(i) * vLoc(j);
        end
    end

    Qtxb = zeros(nx, ny);
    Qtyb = zeros(nx, ny);

end
