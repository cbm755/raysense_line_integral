theta_unscaled = [1 2];
theta = theta_unscaled ./ norm(theta_unscaled)

for N = 2.^[7:18]

  K = 2560;
avg_over = 100;
ss = linspace(0, 1, K);
a = 4; b = 1;  % ref point * 20
Lx = a/20 + ss * theta(1);
Ly = b/20 + ss * theta(2);

g = @(x,y) x.^2 + y.^2;

%% symbolic computation of exact integral
pkg load symbolic
syms t
thetas = sym(theta_unscaled);
thetas = thetas ./ norm(thetas);
lx = sym(a)/20 + t*thetas(1);
ly = sym(b)/20 + t*thetas(2);
exact = int(g(lx,ly), t, [0 1]);
exact = double(exact);

lw = 'linewidth';
plots = false;

err1 = [];
err2 = [];


for ii=1:avg_over
  x = rand(1, N);
  y = rand(1, N);

  if plots
    figure(1); clf;
    plot(x,y,'k.')
    hold on;
    axis equal; axis tight
    plot(Lx, Ly, 'k:', lw, 0.5)
  end

  cloud = [x; y];
  S1 = 0;
  S2 = 0;
  deltar = ss(2) - ss(1);
  for k=1:length(Lx)
    X = [Lx(k); Ly(k)];
    [P, d] = myknn(X, cloud);
    if plots
      plot([X(1) P(1)], [X(2) P(2)], 'r-', lw, 2)
    end
    %% trapezoidal rule weights
    if k == 1 || k == length(Lx)
      w = 1/2;
    else
      w = 1;
    end
    S1 += w*g(X(1), X(2));
    S2 += w*g(P(1), P(2));
  end
  S1 = S1 * deltar;
  S2 = S2 * deltar;

  err1(ii) = abs(S1 - exact);
  err2(ii) = abs(S2 - exact);
end

%20
%err1 = 0.00010958
%err2 = -0.0059666

disp([N K mean(err1) mean(err2) std(err2)])
end
