function [lstar, ustar] = sos_cp(P)
    sdpvar u1 v1 lambda;
    d = 5;
    [s1, s1c] = polynomial([u1, v1], d);
    [s2, s2c] = polynomial([u1, v1], d);
    [s3, s3c] = polynomial([u1, v1], d);
    gi = [u1; v1; 1-u1-v1];
    X = trimapX(P, gi);
    f = X' * X;
    C1 = sos(f - lambda - [s1, s2, s3] * gi);
    C2 = [sos(s1); sos(s2); sos(s3)];
    [C, obj] = sosmodel([C1; C2], -lambda, [], [s1c; s2c; s3c; lambda]);
    optimize(C, obj, []);
    lstar = value(lambda);
    mu = dual(C(2));
    ustar = mu(2:3)/mu(1);
end

function X = trimapX(P, g)
    X = P(:,1)*g(1)^3 + 3*P(:,2)*g(1)^2*g(2) + 3*P(:,3)*g(1)^2*g(3) + ...
        3*P(:,4)*g(1)*g(2)^2 + 6*P(:,5)*g(1)*g(2)*g(3) + ...
        3*P(:,6)*g(1)*g(3)^2 + P(:,7)*g(2)^3 + 3*P(:,8)*g(2)^2*g(3) + ...
        3*P(:,9)*g(2)*g(3)^2 + P(:,10)*g(3)^3;
end