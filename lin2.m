theta_unscaled = [2 3];
theta = theta_unscaled ./ norm(theta_unscaled)
refpt = [1 1]; % ref point * 10

g = @(x,y) cos(y.*x);
%g = @(x,y) x.^2 + y.^2

lw = 'linewidth';
plots = false

%% symbolic computation of exact integral
if is_octave
  pkg load symbolic
  % for knnsearch
  pkg load nan
end
syms t
thetas = sym(theta_unscaled);
thetas = thetas ./ norm(thetas);
lx = sym(refpt(1))/10 + t*thetas(1);
ly = sym(refpt(2))/10 + t*thetas(2);
g(lx,ly)
exact = int(g(lx,ly), t, [0 1])
exact = double(exact)

Nset = 2.^[7:24]
Kset = 10*2.^[0:8]
for Nidx = 1:length(Nset)
  N = Nset(Nidx);
  for Kidx = 1:length(Kset)
    K = Kset(Kidx);
    avg_over = 100;
    ss = linspace(0, 1, K);

    Lx = refpt(1)/10 + ss * theta(1);
    Ly = refpt(2)/10 + ss * theta(2);

    %% trapezoidal rule weights
    trw = ones(1, K);
    trw(1) = 1/2;
    trw(K) = 1/2;

    err1 = [];
    err2 = [];
    tic
    for ii=1:avg_over
      x = rand(1, N);
      y = rand(1, N);

      if plots
	figure(1); clf;
	plot(x, y, 'k.')
	hold on;
	axis equal; axis tight
	plot(Lx, Ly, 'k:', lw, 0.5)
      end

      cloud = [x' y'];
      deltar = ss(2) - ss(1);
      %[P, d] = myknn3(X, cloud);
      X = [Lx' Ly'];
      [idx, d] = knnsearch(cloud, X);
      P = cloud(idx, :);
      if plots
	for k=1:K
	  X = [Lx(k); Ly(k)];
	  Px = [P(k, 1), P(k, 2)];
	  plot([X(1) Px(1)], [X(2) Px(2)], 'r-', lw, 2)
	end
      end
      %% trapezoidal rule
      S1 = sum(trw .* g(Lx, Ly)) * deltar;
      S2 = sum(trw' .* g(P(:,1), P(:,2))) * deltar;

      err1(ii) = abs(S1 - exact);
      err2(ii) = abs(S2 - exact);

      if plots
	disp('paused')
	pause
      end
    end
    tim = toc;
    disp([N K mean(err1) mean(err2) std(err2) tim])
    data{Kidx}(Nidx, :) = [N K mean(err1) mean(err2) std(err2) tim];
  end
  save(sprintf('data_d2new_upto_N%08d', N), 'data')
end
