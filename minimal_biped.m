function minimal_biped

g = 9.81;
m = 50.0;
L_max = 0.5;
D = 0.5;
V = 0.5;
d = D * L_max;
v = V * sqrt(g * L_max);

F_max = 1000;

T = 50;

t_step = d / v;

time = linspace(0, t_step, T);
h = diff(time);

xk = linspace(-0.5 * d, 0.5 * d, T);
yk = 0.5 * ones(1, T);
wk = 0.1 * ones(1, T);
zk = 0.01 * ones(1, T);
Fk = 500 * ones(1, T);

cvx_startup;

cvx_begin;
    
    variable x(T);
    variable y(T);
    variable w(T);
    variable z(T);
    variable F(T);

    minimize sum(F);
    subject to
        x(end) - x(1) == d;          % must move forward by the step length.
        y(1) == y(end);              % periodic motion.
        w(1) == w(end);              % periodic motion.
        z(1) == z(end);              % periodic motion.
        0 <= F <= F_max;             % leg strength.
        x.^2 + y.^2 <= L_max^2;      % leg length.
        for i = 1:(T-1)
            Lk = sqrt(xk(i+1)^2 + yk(i+1)^2);

            % affine approximation for x dynamics.
            cxk = wk(i+1) - wk(i) - h(i) * ( Fk(i+1) / m * xk(i+1) / Lk );
            % w(i+1)  w(i)  x(i+1)  y(i+1)  F(i+1)
            Gxk = [1, ...
                   -1, ...
                   -h(i) * Fk(i+1) / m * ( -xk(i+1)^2 / Lk^3 + 1 / Lk ), ...
                   -h(i) * Fk(i+1) / m * xk(i+1) * ( -yk(i+1) / Lk^3 ), ...
                   -h(i) * xk(i+1) / (m * Lk)];

            % affine approximation for y dynamics.
            cyk = zk(i+1) - zk(i) - h(i) * ( Fk(i+1) / m * yk(i+1) / Lk - g );
            % z(i+1)  z(i)  x(i+1)  y(i+1)  F(i+1)
            Gyk = [1, ...
                   -1, ...
                   -h(i) * Fk(i+1) / m * yk(i+1) * ( -xk(i+1) / Lk^3 ), ...
                   -h(i) * Fk(i+1) / m * ( - yk(i+1)^2 / Lk^3 + 1 / Lk ), ...
                   -h(i) * yk(i+1) / (m * Lk)];

            % non-convex dynamics constraints.
            xvec = [w(i+1); w(i); x(i+1); y(i+1); F(i+1)];
            xveck = [wk(i+1); wk(i); xk(i+1); yk(i+1); Fk(i+1)];
            cxk + Gxk * (xvec - xveck)  == 0;

            yvec = [z(i+1); z(i); x(i+1); y(i+1); F(i+1)];
            yveck = [zk(i+1); zk(i); xk(i+1); yk(i+1); Fk(i+1)];
            cyk + Gyk * (yvec - yveck) == 0;

            % convex dynamics constraints.
            x(i+1) - x(i) - h(i) * w(i+1) == 0;
            y(i+1) - y(i) - h(i) * z(i+1) == 0;

        end

cvx_end;

subplot(1, 2, 1);
plot(time, F);
xlabel('time (s)');
ylabel('force (N)');
subplot(1, 2, 2);
L = sqrt(x.^2 + y.^2);
plot(time, L);
xlabel('time (s)');
ylabel('length (m)');
print -dsvg minimal_biped_cvx.svg

end
