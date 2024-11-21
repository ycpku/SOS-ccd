function [t1, t2, v1, v2] = generate_bezier_patch(n)
    num = (n+1)^2;
    t1 = randn(3,num);
    t2 = randn(3,num);
    v1 = randn(3,num);
    v2 = randn(3,num);
    shift = [0;0;1/sqrt(2)];
    t1 = t1 + shift;
    t2 = t2 - shift;
    v1 = v1 - shift;
    v2 = v2 + shift;
end