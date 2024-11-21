%% without culling
num_tests = 1000;

for order = [1 2 3]
for type = ["triangle" "patch"]
tic
if(strcmp(type,"triangle"))
    num_cp = (order+1)*(order+2)/2;
elseif(strcmp(type,"patch"))
    num_cp = (order+1)*(order+1);
else
    assert(false)
end

data(1:num_tests) = struct('p1',zeros(3,num_cp),'p2',zeros(3,num_cp), ...
    'v1',zeros(3,num_cp),'v2',zeros(3,num_cp),'t',[],'u',zeros(1,5), ...
    'yalmiptime',[],'solvertime',[],'info',[],'problem',[],'exact_recovery',0);
for i = 1:num_tests
    % [P1, P2, V1, V2] = generate_bezier_triangle(order);
    [P1, P2, V1, V2] = read_from_file("data/"+type+int2str(order)+"/case"+int2str(i)+".dat",type,order);
    data(i).p1 = P1;data(i).p2 = P2;data(i).v1 = V1;data(i).v2 = V2;
    %[t, u, diagnostics, eig_vals] = ccd_triangle(P1, P2, V1, V2, order);
    if(strcmp(type,"triangle"))
        [t, u, diagnostics, eig_vals] = ccd_triangle_speedup(P1, P2, V1, V2, order);
    elseif(strcmp(type,"patch"))
        [t, u, diagnostics, eig_vals] = ccd_patch(P1, P2, V1, V2, order);
    else
        assert(false)
    end
    data(i).t = t;
    data(i).u = u;
    data(i).yalmiptime = diagnostics.yalmiptime;
    data(i).solvertime = diagnostics.solvertime;
    data(i).info = diagnostics.info;
    data(i).problem = diagnostics.problem;
    if eig_vals(end-1)<1e-4
        data(i).exact_recovery = 1;
    end
end
timeElapsed = toc
save(type+int2str(order)+".mat","data","timeElapsed")
end
end
%% culling
num_tests = 10;
order = 3;
eps_lambda = 1e-3;
eps_gamma = 1e-3;
eps_d = 1e-6;
type = 'triangle';

if(strcmp(type,'triangle'))
    num_cp = (order+1)*(order+2)/2;
elseif(strcmp(type,'patch'))
    num_cp = (order+1)*(order+1);
else
    assert(false)
end

data(1:num_tests) = struct('p1',zeros(3,num_cp),'p2',zeros(3,num_cp), ...
    'v1',zeros(3,num_cp),'v2',zeros(3,num_cp),'t',[],'u',zeros(1,5), ...
    'yalmiptime',[],'solvertime',[],'info',[],'problem',[],'exact_recovery',0);
for i = 1:num_tests
    [P1, P2, V1, V2] = generate_bezier_triangle(order);
    data(i).p1 = P1;data(i).p2 = P2;data(i).v1 = V1;data(i).v2 = V2;
    %[t, u, diagnostics, eig_vals] = ccd_triangle(P1, P2, V1, V2, order);
    if(strcmp(type,'triangle'))
        [~, u, diagnostics, eig_vals] = ip_triangle(P1, P2, order);
        if eig_vals(end-1)<eps_lambda
            data(i).t = 0;
            data(i).u = u;
            data(i).yalmiptime = diagnostics.yalmiptime;
            data(i).sovlertime = diagnostics.solvertime;
            data(i).info = diagnostics.info;
            data(i).problem = diagnostics.problem;
            continue;
        end
        [~, u, diagnostics, eig_vals] = ip_triangle(P2, P1, order);
        if eig_vals(end-1)<eps_lambda
            data(i).t = 0;
            data(i).u = u;
            data(i).yalmiptime = diagnostics.yalmiptime;
            data(i).sovlertime = diagnostics.solvertime;
            data(i).info = diagnostics.info;
            data(i).problem = diagnostics.problem;
            continue;
        end
        [gamma, u, diagnostics, eig_vals] = nc_triangle(P1, P2, V1, V2, order);
        if gamma>=eps_gamma
            data(i).t = -1;
            data(i).u = u;
            data(i).yalmiptime = diagnostics.yalmiptime;
            data(i).sovlertime = diagnostics.solvertime;
            data(i).info = diagnostics.info;
            data(i).problem = diagnostics.problem;
            continue;
        end
        [t, u, diagnostics, eig_vals] = ccd_triangle_speedup(P1, P2, V1, V2, order);
    elseif(strcmp(type,'patch'))
        [t, u, diagnostics, eig_vals] = ccd_patch(P1, P2, V1, V2, order);
    else
        assert(false)
    end
    data(i).t = t;
    data(i).u = u;
    data(i).yalmiptime = diagnostics.yalmiptime;
    data(i).sovlertime = diagnostics.solvertime;
    data(i).info = diagnostics.info;
    data(i).problem = diagnostics.problem;
    if eig_vals(end-1)<1e-4
        data(i).exact_recovery = 1;
    end
end
