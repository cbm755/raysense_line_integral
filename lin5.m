theta_unscaled = [2 3 4 5 6];
theta = theta_unscaled ./ norm(theta_unscaled)
refpt = [1 1 1 1 1]; % ref point * 10

g = @(x,y,z,w,w2) cos(y.*x) - sin(z).*w.*w2;

%% symbolic computation of exact integral
if is_octave
  pkg load symbolic
  % for knnsearch
  pkg load nan
end
syms t
thetas = sym(theta_unscaled);
thetas = thetas ./ norm(thetas);
d = 5
for i=1:d
  l{i} = sym(refpt(i))/10 + t*thetas(i);
end
g(l{1},l{2},l{3},l{4},l{5})
exact = int(g(l{1},l{2},l{3},l{4},l{5}), t, [0 1])
exact = double(exact)

Nset = 2.^[7:23]
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
    Lw = refpt(4)/10 + ss * theta(4);
    Lw2 = refpt(5)/10 + ss * theta(5);

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
      z = rand(1, N);
      w = rand(1, N);
      w2 = rand(1, N);

      cloud = [x' y' z' w' w2'];
      deltar = ss(2) - ss(1);
      %[P, d] = myknn3(X, cloud);
      X = [Lx' Ly' Lz' Lw' Lw2'];
      [idx, d] = knnsearch(cloud, X);
      P = cloud(idx, :);
      %% trapezoidal rule
      S1 = sum(trw .* g(Lx, Ly, Lz, Lw, Lw2)) * deltar;
      S2 = sum(trw' .* g(P(:,1), P(:,2), P(:,3), P(:,4), P(:,5))) * deltar;

      err1(ii) = abs(S1 - exact);
      err2(ii) = abs(S2 - exact);
    end
    tim = toc;
    disp([N K mean(err1) mean(err2) std(err2) tim])
    data{Kidx}(Nidx, :) = [N K mean(err1) mean(err2) std(err2) tim];
  end
  save(sprintf('data_d5new_upto_N%08d', N), 'data')
end
