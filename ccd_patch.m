 function [lstar, ustar, diagnostics, eig_vals] = ccd_patch(P1, P2, V1, V2, n)
    sdpvar u1 v1 u2 v2 t lambda;
    if(n==1)
        d1 = 4;
        d2 = 2; % mixed degree
    elseif(n==2)
        d1 = 4;
        d2 = 3; % mixed degree
    else
        d1 = 6;
        d2 = 4; % mixed degree
    end
    [s1, s1c] = polynomial([u1, v1, u2, v2, t], d1);
    [s2, s2c] = polynomial([u1, v1, u2, v2, t], d1);
    [s3, s3c] = polynomial([u1, v1, u2, v2, t], d1);
    [s4, s4c] = polynomial([u1, v1, u2, v2, t], d1);
    [s5, s5c] = polynomial([u1, v1, u2, v2, t], d1);
    [p1, p1c] = polynomial([u1, v1, u2, v2, t], d2);
    [p2, p2c] = polynomial([u1, v1, u2, v2, t], d2);
    [p3, p3c] = polynomial([u1, v1, u2, v2, t], d2);
    gi = [u1*(1-u1); v1*(1-v1); u2*(1-u2); v2*(1-v2); t*(1-t)]; % higher degree description
    X1 = quadmapX(P1, u1, v1, n) + t * quadmapX(V1, u1, v1, n);
    X2 = quadmapX(P2, u2, v2, n) + t * quadmapX(V2, u2, v2, n);
    hi = X1 - X2;
    C1 = sos(t - lambda - [p1, p2, p3] * hi - [s1, s2, s3, s4, s5] * gi);
    C2 = [sos(s1); sos(s2); sos(s3); sos(s4); sos(s5)];
    C3 = [0<=lambda; lambda<=1]; % constraining lambda
    [C, obj] = sosmodel([C1; C2; C3], -lambda, ...
        [], [s1c; s2c; s3c; s4c; s5c; p1c; p2c; p3c; lambda]);
    diagnostics = optimize(C, obj, []);
    lstar = value(lambda);
    mu = dual(C(end-5));
    eig_vals = eig(mu);
    ustar = mu(2:6)/mu(1);
    ustar = ustar(end:-1:1);
end

function X = quadmapX(P, u, v, n)
if(n==1)
    X = P(:,1)*(1-u)*(1-v) + P(:,2)*(1-u)*v + P(:,3)*u*(1-v) + P(:,4)*u*v;
elseif(n==2)
    X = P(:,1)*(1-u)^2*(1-v)^2 + 2*P(:,2)*(1-u)^2*v*(1-v) + P(:,3)*(1-u)^2*v^2 + ...
        2*P(:,4)*u*(1-u)*(1-v)^2 + 4*P(:,5)*u*(1-u)*v*(1-v) + 2*P(:,6)*u*(1-u)*v^2 + ...
        P(:,7)*u^2*(1-v)^2 + 2*P(:,8)*u^2*v*(1-v) + P(:,9)*u^2*v^2;
else
    X = P(:,1)*(1-u)^3*(1-v)^3 + 3*P(:,2)*(1-u)^3*v*(1-v)^2 + 3*P(:,3)*(1-u)^3*v^2*(1-v) + P(:,4)*(1-u)^3*v^3 + ...
        3*P(:,5)*u*(1-u)^2*(1-v)^3 + 9*P(:,6)*u*(1-u)^2*v*(1-v)^2 + 9*P(:,7)*u*(1-u)^2*v^2*(1-v) + 3*P(:,8)*u*(1-u)^2*v^3 + ...
        3*P(:,9)*u^2*(1-u)*(1-v)^3 + 9*P(:,10)*u^2*(1-u)*v*(1-v)^2 + 9*P(:,11)*u^2*(1-u)*v^2*(1-v) + 3*P(:,12)*u^2*(1-u)*v^3 + ...
        P(:,13)*u^3*(1-v)^3 + 3*P(:,14)*u^3*v*(1-v)^2 + 3*P(:,15)*u^3*v^2*(1-v) + P(:,16)*u^3*v^3;
end
end