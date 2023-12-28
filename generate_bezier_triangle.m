function [t1, t2, v1, v2] = generate_bezier_triangle()
    t1 = randn(3,10);
    t2 = randn(3,10);
    v1 = randn(3,10);
    v2 = randn(3,10);
    shift = [0;0;1/sqrt(2)];
    t1 = t1 + shift;
    t2 = t2 - shift;
    v1 = v1 - shift;
    v2 = v2 + shift;
end