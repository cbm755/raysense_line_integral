%%
% TODO:
% aspect ratio hacks: instead I just set last three windows to be 55% or 60% there usual width

for d=2:5
  N = 2^24
  N2 = 2^23
  if d == 2
    load(sprintf('data/data_d2new_upto_N%08d.mat', N))
    %load('nonunif_d2_upto_N16776417')
  elseif d == 3
    load(sprintf('data/data_d3new_upto_N%08d.mat', N))
  elseif d == 4
    load(sprintf('data/data_d4new_upto_N%08d.mat', N))
  elseif d == 5
    load(sprintf('data/data_d5new_upto_N%08d.mat', N2))
  else
    error('no data')
  end

  figure(3+d); clf;
  A = data;

lwdot = 1.4
lwdata = 1.5
lwpred = 2.5
lwthm = 1.4

leg = {}
colours = {}
legends_handles = []
labels = {}
lw = 'linewidth';

for i=1:4
  if i == 1
    hquad = loglog(A{i}(:,1), A{i}(:,3), 'k--', lw, 0.7);
    hold on
  else
    plot(A{i}(:,1), A{i}(:,3), 'k--', lw, 0.7);
  end
  K = A{i}(1,2);
  %set(h, 'color', colours{i})
  %labels{end+1} = ['k = ' num2str(K) ' quad'];
end

for i=1:length(A)
  h = plot(A{i}(:,1), A{i}(:,4), lw, lwdata);
  K = A{i}(1,2);
  legends_handles(i) = h
  labels{i} = ['n_r = ' num2str(K)];
  colours{i} = get(h, 'color');
end
set(gca, 'fontsize', 14)

%% we only label one of the quadrature ones
legends_handles(end+1) = hquad;
labels{end+1} = ['quadrature'];

N = A{i}(:,1);
M = N.^(1/d)
xlabel('N')
if d == 2
  h = plot(N, 0.02*N.^(-1), 'k:', lw, lwdot)
  legends_handles(end+1) = h;
  labels{end+1} = 'N^{-1}';
  h = plot(N, 0.0025*N.^(-0.5), 'k:', lw, lwdot)
  legends_handles(end+1) = h;
  labels{end+1} = 'N^{-0.5}';
  h = plot(N, 0.0052*(log(N)./N), 'k--', lw, lwpred)
  legends_handles(end+1) = h;
  labels{end+1} = 'ln(N)/N';
  %plot(N, 0.002*M.*(log(N)./N), 'k-.', lw, 2)
  %plot(N, 0.002*(log(N)./sqrt(N)), 'k-.', lw, 2)
  h = plot(N, 0.02*(log(N).^(1/2)./sqrt(N)), 'k-.', lw, lwthm)
  legends_handles(end+1) = h;
  labels{end+1} = '(ln(N)/N)^{1/2}';
  %plot(N, 0.033*N.^(-0.6), 'r:', lw, 3)
  %legends_handles(end+1) = h1;x  %labels{end+1} = 'N^{-1/2}';
  %labels{end+1} = 'ln(N)/(N^{1/2})';

  legend(legends_handles, labels, 'location', 'eastoutside')
  ylim([1e-7 4e-3])
  xlim([1e2 2.1e6])
  ylabel('error')
elseif d == 3
  plot(N, 0.05*N.^(-0.5), 'k:', lw, lwdot)
  plot(N, 0.4*N.^(-1), 'k:', lw, lwdot)
  h = plot(N, 0.033*(log(N)./N).^(2/3), 'k--', lw, lwpred)
  %h2 = plot(N, 0.03*M.*(log(N)./N).^(2/3), 'k-.', lw, 2)
  h2 = plot(N, 0.055*(log(N)./N).^(1/3), 'k-.', lw, lwthm)
  %h3 = plot(N, 0.08*N.^(-0.5), 'r:-', lw, 2)
  %labels{end+1} = 'N^{-1/2}';
  %labels{end+1} = '(log(N)/N)^{2/3}';
  %labels{end+1} = 'N^{-1}';
  %legend(labels, 'location', 'eastoutside')
  legs = {'(ln(N)/N)^{2/3}', '(ln(N)/N)^{1/3}'};
  legs{end+1} = 'N^{-0.5}';
  legend([h h2], legs);
  ylim([1e-6 2e-2])
  xlim([1e2 3e7])
elseif d == 4
  plot(N, 0.45*N.^(-1), 'k:', lw, lwdot)
  plot(N, 0.04*N.^(-0.5), 'k:', lw, lwdot)
  h = plot(N, 0.02*(log(N)./N).^(2/4), 'k--', lw, lwpred)
  h2 = plot(N, 0.036*(log(N)./N).^(1/4), 'k-.', lw, lwthm)
  %h3 = plot(N, 0.01*(log(N).^(1/4))./(N.^(1/4)), 'k-.', lw, 2)
  %h3 = plot(N, 0.04*N.^(-0.4), 'r:-', lw, 2)
  %labels{end+1} = 'N^{-1/2}';
  %labels{end+1} = '(log(N)/N)^{1/2}';
  %labels{end+1} = 'N^{-1}';
  legs = {'(ln(N)/N)^{1/2}', '(ln(N)/N)^{1/4}'};
  legs{end+1} = 'N^{-0.4}';
  legend([h h2], legs)
  ylim([1e-5 2e-2])
  xlim([1e2 3e7])
elseif d == 5
  plot(N, 0.3*N.^(-1), 'k:', lw, lwdot)
  plot(N, 0.028*N.^(-0.5), 'k:', lw, lwdot)
  h = plot(N, 0.01*(log(N)./N).^(2/5), 'k--', lw, lwpred)
  h2 = plot(N, 0.024*(log(N)./N).^(1/5), 'k-.', lw, lwthm)
  %h3 = plot(N, 0.01*N.^(-0.3), 'r:-', lw, 2)
  %labels{end+1} = 'N^{-1/2}';
  %labels{end+1} = '(log(N)/N)^{1/2}';
  %labels{end+1} = 'N^{-1}';
  legs = {'(ln(N)/N)^{2/5}', '(ln(N)/N)^{1/5}'};
  %legs{end+1} = 'N^{-0.3}';
  legend([h h2], legs);
  ylim([1e-5 2e-2])
  xlim([1e2 1e7])
else
  error('invalid dimension')
end
title(['d = ' num2str(d)])
%xlim([100 300000])
b = ['line_integral_uniform_d' num2str(d)']
print(b, '-dpng')
print(b, '-depsc2')
assert(~ system(['epstopdf ' b '.eps']))
end
