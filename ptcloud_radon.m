function [R1, R2] = ptcloud_radon(cloud, nr, theta_set, Mbins, g, plots)
%PTCLOUD_RADON   Approx radon transform on a point cloud
%  [R1, R2] = ptcld_radon(CLOUD, NR, THETA_SET, MBINS, G)
%     CLOUD: 2 x N: a 2-d point cloud
%     NR: samples per unit-length ray
%     THETA_SET: defaults to 0:1:180
%     MBINS: number of bins in transform
%     G: function handle g(x,y)
%        This is the function data to be transformed.
%        TODO: instead it should be taken as
%        attributes of CLOUD or another vector.
%
%  An extra boolean argument plots is a hack for making debugging plots.


[d, N] = size(cloud);

assert(d == 2, 'only works in 2-d for now')

if nargin < 3
  theta_set = 0:1:180;
end
if nargin < 4
  Mbins = 31;
end
if nargin < 6
  plots = false;
end
Nth = length(theta_set);

%% config, stuff for debugging
lw = 'linewidth';
sym_exact = false;

%% point cloud
x = cloud(1, :);
y = cloud(2, :);

if plots
  figure(1); clf;
  plot(x,y,'k.', lw, 0.5, 'markersize', 5)
  hold on;
  axis equal; axis tight
end

R1 = zeros(Mbins, Nth);
R2 = zeros(Mbins, Nth);

Time = 0
for theta_idx = 1:length(theta_set)
  tic
  theta = theta_set(theta_idx)
  perp = [cosd(theta); sind(theta)];
  dirvec = [-sind(theta); cosd(theta)];

  if sym_exact
    pkg load symbolic
    dirvecs = [cosd(sym(theta)); sind(sym(theta))]
    perps = [-sind(sym(theta)); cosd(sym(theta))]
  end

  if plots
    [xx,yy] = meshgrid(linspace(min(x), max(x), 512), linspace(min(y), max(y), 512));
    figure(2); clf;
    pcolor(xx,yy,g(xx,yy))
    shading flat
    hold on
    plot(x,y,'k.', lw, 0.5, 'markersize', 4)
    hold on;
    axis equal; axis tight
    set(gca, 'clim', [0 1.2])
    colormap(flipud(pink))
  end

  if plots
    offsetx = dirvec(1) + ((0:(Mbins-1))./(Mbins-1) - 0.5)*perp(1);
    offsety = dirvec(2) + ((0:(Mbins-1))./(Mbins-1) - 0.5)*perp(2);
    plot(offsetx, offsety, 'b-', lw, 2)
  end

  for idx=1:Mbins
    ps = (idx-1)/(Mbins-1);
    offset = [0; 0] + (ps - 0.5)*perp;
    ss = linspace(-0.5, 0.5, nr);
    Lx = offset(1) + ss * dirvec(1);
    Ly = offset(2) + ss * dirvec(2);
    if plots && theta_idx == 5
      plot(Lx, Ly, 'b:', lw, 1.5)
    end

    %% optional symbolic computation of exact integral
    if sym_exact
      pss = sym(idx-1)/(Mbins-1);
      offsets = (pss - sym(1)/2)*perps;
      pkg load symbolic
      t = sym('t')
      lx = offsets(1) + t*dirvecs(1);
      ly = offsets(2) + t*dirvecs(2);
      exact = int(g(lx,ly), t, [-sym(1)/2 sym(1)/2]);
      exact = double(exact);
    end

    %% quadrature
    S1 = 0;
    S2 = 0;
    deltar = ss(2) - ss(1);
    Px = zeros(size(Lx));
    Py = zeros(size(Lx));
    for k=1:length(Lx)
      X = [Lx(k); Ly(k)];
      % TODO: replace with faster implementation
      [P, d] = myknn(X, cloud);
      if plots && theta_idx == 5
	plot(P(1), P(2), 'bo', lw, 1, 'markersize', 4)
	plot([X(1) P(1)], [X(2) P(2)], 'b:', lw, 1, 'markersize', 4)
      end
      Px(k) = P(1);
      Py(k) = P(2);
      %% trapezoidal rule weights
      if k == 1 || k == length(Lx)
	w = 1/2;
      else
	w = 1;
      end
      S1 = S1 + w*g(X(1), X(2));
      S2 = S2 + w*g(P(1), P(2));
    end
    if plots && theta_idx == 5
      plot(Px, Py, 'bo', lw, 1.5)
    end
    S1 = S1 * deltar;
    S2 = S2 * deltar;

    R1(idx, theta_idx) = S1;
    R2(idx, theta_idx) = S2;

    if sym_exact
      err1(end+1) = abs(S1 - exact);
      err2(end+1) = abs(S2 - exact);
    end
    if plots
      %drawnow
    end
  end
  if plots && theta_idx == 5
    pause
    print(['example_cp_M' num2str(Mbins) '.png'], '-dpng')
  end
  Time += toc
end
