num_tests = 10;
data(1:num_tests) = struct('p1',zeros(3,16),'p2',zeros(3,16), ...
    'v1',zeros(3,16),'v2',zeros(3,16),'t',[],'u',zeros(1,5), ...
    'yalmiptime',[],'sovlertime',[],'info',[],'problem',[],'exact_recovery',0);
for i = 1:num_tests
    [P1, P2, V1, V2] = generate_bezier_triangle();
    data(i).p1 = P1;data(i).p2 = P2;data(i).v1 = V1;data(i).v2 = V2;
    [t, u, diagnostics, eig_vals] = sos_ccd(P1, P2, V1, V2);
    data(i).t = t;
    data(i).u = u;
    data(i).yalmiptime = diagnostics.yalmiptime;
    data(i).sovlertime = diagnostics.solvertime;
    data(i).info = diagnostics.info;
    data(i).problem = diagnostics.problem;
    if eig_vals(end-1)<1e-6
        data(i).exact_recovery = 1;
    end
end