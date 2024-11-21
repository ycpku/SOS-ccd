function [lstar, ustar, diagnostics, eig_vals] = ip_triangle(P1, P2, n)
    sdpvar u1 v1 u2 v2 t lambda;
    if(n==1)
        d1 = 2;
        d2 = 2;
    elseif(n==2)
        d1 = 2;
        d2 = 2;
    else
        d1 = 4;
        d2 = 2; % mixed degree
    end
    [s1, s1c] = polynomial([u1, v1, u2, v2, t], d1);
    [s2, s2c] = polynomial([u1, v1, u2, v2, t], d1);
    [s3, s3c] = polynomial([u1, v1, u2, v2, t], d1);
    [s4, s4c] = polynomial([u1, v1, u2, v2, t], d1);
    [s5, s5c] = polynomial([u1, v1, u2, v2, t], d1);
    [s6, s6c] = polynomial([u1, v1, u2, v2, t], d1);
    [p1, p1c] = polynomial([u1, v1, u2, v2, t], d2);
    [p2, p2c] = polynomial([u1, v1, u2, v2, t], d2);
    [p3, p3c] = polynomial([u1, v1, u2, v2, t], d2);
    gi = [u1; v1; 1-u1-v1; u2; v2; 1-u2-v2];
    X1 = trimapX(P1, gi(1:3),n);
    X2 = trimapX(P2, gi(4:6),n);
    hi = X1 - X2;
    C1 = sos(X1(1) - lambda - [p1, p2, p3] * hi - [s1, s2, s3, s4, s5, s6] * gi);
    C2 = [sos(s1); sos(s2); sos(s3); sos(s4); sos(s5); sos(s6)];
    [C, obj] = sosmodel([C1; C2], -lambda, ...
        [], [s1c; s2c; s3c; s4c; s5c; s6c; p1c; p2c; p3c; lambda]);
    diagnostics = optimize(C, obj, []);
    lstar = value(lambda);
    mu = dual(C(2));
    eig_vals = eig(mu);
    ustar = mu(2:5)/mu(1);
    ustar = ustar(end:-1:1);
end

function X = trimapX(P, g, n)
if(n==1)
    X = P(:,1)*g(1) + P(:,2)*g(2) + P(:,3)*g(3);
elseif(n==2)
    X = P(:,1)*g(1)^2 + 2*P(:,2)*g(1)*g(2) + 2*P(:,3)*g(1)*g(3) + ...
        P(:,4)*g(2)^2 + 2*P(:,5)*g(2)*g(3) + P(:,6)*g(3)^2;
else
    X = P(:,1)*g(1)^3 + 3*P(:,2)*g(1)^2*g(2) + 3*P(:,3)*g(1)^2*g(3) + ...
        3*P(:,4)*g(1)*g(2)^2 + 6*P(:,5)*g(1)*g(2)*g(3) + ...
        3*P(:,6)*g(1)*g(3)^2 + P(:,7)*g(2)^3 + 3*P(:,8)*g(2)^2*g(3) + ...
        3*P(:,9)*g(2)*g(3)^2 + P(:,10)*g(3)^3;
end
end