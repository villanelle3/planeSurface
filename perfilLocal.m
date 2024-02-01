function [uLoc, vLoc] = perfilLocal(xi, ny, y, umax, ni)
    eta = zeros(1, ny);
    g = zeros(1, ny);
    f = zeros(1, ny);

    if (xi <= 0)
        uLoc = umax;
        vLoc = zeros(size(umax));
    else
        for j = 1:ny
            eta(j) = y(j) * sqrt(umax / (ni * xi));

            if (eta(j) >= 0 && eta(j) <= 5)
                g(j) = (0.332057 * eta(j) + 0.00059069 * eta(j)^4 + 0.00000288 * ...
                    eta(j)^5 * exp(eta(j)^2 / 4 - 1)) / (1 + 0.00869674 * eta(j)^3 ...
                    + 0.00000288 * eta(j)^5 * exp(eta(j)^2 / 4 - 1));
                f(j) = 1 / 6.0 * eta(j)^2.0 * tanh((4 * pi)^(1 / 3) * ...
                    (atan(sqrt(4 * pi / (eta(j)^3))))^(2 / 3));
            elseif (eta(j) > 5)
                g(j) = (0.332057 * eta(j) + 0.00059069 * eta(j)^4 + 0.00000288 * ...
                    eta(j)^5 * exp(eta(j)^2 / 4 - 1)) / (1 + 0.00869674 * eta(j)^3 ...
                    + 0.00000288 * eta(j)^5 * exp(eta(j)^2 / 4 - 1));
                f(j) = eta(j) - 1.7208;
            end
        end

        uLoc = zeros(1, ny);
        vLoc = zeros(1, ny);

        for j = 1:ny
            if (eta(j) > 0 && eta(j) < 5)
                vLoc(j) = 0.5 * sqrt(ni * umax / xi) * (eta(j) * g(j) - f(j));
                uLoc(j) = g(j) * umax;
            elseif (eta(j) > 5)
                vLoc(j) = 0.5 * sqrt(ni * umax / xi) * (eta(j) * g(j) - f(j));
                uLoc(j) = g(j) * umax;
            end
        end
    end
end
