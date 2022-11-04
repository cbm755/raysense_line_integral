theta_unscaled = [2 3 4];
theta = theta_unscaled ./ norm(theta_unscaled)
refpt = [1 1 1]; % ref point * 10

g = @(x,y,z) x.^3 + cos(y) - sin(z);
g = @(x,y,z) cos(y.*x) - sin(z);

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
lz = sym(refpt(3))/10 + t*thetas(3);
g(lx,ly,lz)
exact = int(g(lx,ly,lz), t, [0 1])
exact = double(exact)

Nset = 2.^[7:24]
Kset = 10*2.^[0:8]
for Nidx = 1:length(Nset)
  N = Nset(Nidx);
  for Kidx = 1:length(Kset)
    K = Kset(Kidx);
    avg_over = 50;
    ss = linspace(0, 1, K);

    Lx = refpt(1)/10 + ss * theta(1);
    Ly = refpt(2)/10 + ss * theta(2);
    Lz = refpt(3)/10 + ss * theta(3);

    err1 = [];
    err2 = [];
    tic
    for ii=1:avg_over
      x = rand(1, N);
      y = rand(1, N);
      z = rand(1, N);

      if plots
	figure(1); clf;
	plot3(x,y,z,'k.')
	hold on;
	axis equal; axis tight
	plot3(Lx, Ly, Lz, 'k:', lw, 0.5)
      end

      cloud = [x' y' z'];
      deltar = ss(2) - ss(1);
      %[P, d] = myknn3(X, cloud);
      X = [Lx' Ly' Lz'];
      [idx, d] = knnsearch(cloud, X);
      P = cloud(idx, :);
      if plots
	for k=1:length(Lx)
	  X = [Lx(k); Ly(k); Lz(k)];
	  Px = [P(k, 1), P(k, 2), P(k, 3)];
	  plot3([X(1) Px(1)], [X(2) Px(2)], [X(3) Px(3)], 'r-', lw, 2)
	end
      end
      %% trapezoidal rule weights
      w = ones(1, K);
      w(1) = 1/2;
      w(K) = 1/2;
      S1 = sum(w .* g(Lx, Ly, Lz)) * deltar;
      S2 = sum(w' .* g(P(:,1), P(:,2), P(:,3))) * deltar;

      err1(ii) = abs(S1 - exact);
      err2(ii) = abs(S2 - exact);
    end
    tim = toc;
    disp([N K mean(err1) mean(err2) std(err2) tim])
    data{Kidx}(Nidx, :) = [N K mean(err1) mean(err2) std(err2) tim];
  end
  save(sprintf('data_d3new_upto_N%08d', N), 'data')
end
