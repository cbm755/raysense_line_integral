%%
% Tomography experiment using approximate point cloud radon transform

%% Octave packages
% for iradon (filtered back projection)
if is_octave()
  pkg load image
  % for optional symbolic debugging
  %pkg load symbolic
end

% if false, use the existing cloud, else make new one
new_cloud = true

if new_cloud
  N = 30000
  K = 64;
  %% point cloud
  x = 1.44*rand(1, N) - 0.72;
  y = 1.44*rand(1, N) - 0.72;

  % use a third dimension to threshold against a density function
  z = rand(1, N);
  % rho1 = @(x,y) (sin(2*pi*x).*sin(2*pi*y) + 1.25)/2.25;
  % rho = @(x,y) rho1(x,y) .* (abs(x) <= 0.5) .* (abs(y) <= 0.5) + ...
  %       0.5*(abs(y) > 0.5 | abs(x) > 0.5);
  % 				%0.5*(abs(x) > 0.5) .* ...
  rho = @(x,y) -3*x.*exp(-9*(x.^2 +y.^2)) + 0.5
  I = z < rho(x, y);
  x = x(I);
  y = y(I);

  cloud = [x; y];
  N = length(x)

  uniq = randi(99999);
end


dt = 1;
theta_set = 0:dt:180;
Nth = length(theta_set);
Mbins = 101

%% config
lw = 'linewidth';
plots = true;


if plots
  figure(1); clf;
  plot(x,y,'k.', lw, 0.5, 'markersize', 5)
  hold on;
  axis equal; axis tight
  print(sprintf('point_cloud_%05d.png', uniq), '-dpng')
end

%% density
g = @(x,y) 1.0*(...
		 ((((x + 0.22).^2 + (y - 0.12).^2) <= 0.15^2) | ...
 		  (((x - 0.22).^2 + (y - 0.12).^2) <= 0.12^2) | ...
		  (((x - 0).^2 + (y + 0.15).^2) <= 0.2^2)) & ...
		 ((x.^2 + y.^2) >= 0.25^2) ...
	       ) + ...
    0.5*((x.^2 + y.^2) <= 0.2^2) ...
;
[xx,yy] = meshgrid(linspace(-0.72, 0.72, 512), linspace(-0.72, 0.72, 512));

base = sprintf('%05d_M%d_K%d_dt%d', uniq, Mbins, K, dt);

if plots
  figure(10); clf;
  pcolor(xx, yy, rho(xx, yy))
  axis tight; axis equal
  shading flat
  %colormap(flipud(viridis))
  colormap(flipud(gray))
  print(['ptcloud_density_' base '.png'], '-dpng')
end

disp('paused')
pause

[R1, R2] = ptcloud_radon(cloud, K, 0:15:180, 21, g, true);

[R1, R2] = ptcloud_radon(cloud, K, theta_set, Mbins, g);


save(['sav_' base], 'cloud', 'R1', 'R2', 'N', 'theta_set', 'Mbins', 'K', 'g')

figure(4); clf;
pcolor(R1)
xlabel('theta')
title('quadrature')
figure(5); clf;
pcolor(R2)
xlabel('\theta')
%title('raysense quadrature')
shading flat
colormap(flipud(hot))
print(['sonogram_' base '.png'], '-dpng')


pad = 10;
R1p = [zeros(pad,Nth); R1; zeros(pad,Nth)];
R2p = [zeros(pad,Nth); R2; zeros(pad,Nth)];

figure(8); clf;
A = iradon(R1p);
A = flipud(A);
pcolor(A)
title('filtered back projection: quadrature, padded')
axis equal; axis tight

figure(9); clf;
A = iradon(R2p, theta_set, 'linear', 'Hamming'); % 'Shepp-Logan');
A = flipud(A);
pcolor(A)
%title('filtered back projection: RaySense, padded')
axis equal; axis tight
shading flat
colormap(flipud(pink))
print(['fbp_' base '.png'], '-dpng')

figure(10); clf;
A = iradon(R2p, theta_set, 'linear', 'none');
A = flipud(A);
pcolor(A)
title('filtered back projection: RaySense, padded')
axis equal; axis tight
shading flat
colormap(flipud(pink))
