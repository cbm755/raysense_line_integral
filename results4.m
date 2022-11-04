for d=2:4
  N = 2^21
  N2 = 2^21
  if d == 2
    load(sprintf('data_d2new_upto_N%08d.mat', N))
  elseif d == 3
    %load(sprintf('data_upto_N%08d.mat', N))
    load(sprintf('data_d3new_upto_N%08d.mat', N))
  elseif d == 4
    %load(sprintf('data_d4_upto_N%08d.mat', N))
    load(sprintf('data_d4new_upto_N%08d.mat', N2))
  else
    error('no data')
  end

  figure(3+d); clf;
  A = data;

leg = {}
colours = {}
labels = {}
lw = 'linewidth';
for i=1:length(A)
  if i == 1
    h = loglog(A{i}(:,1), A{i}(:,4), lw, 1.5);
    hold on
  else
    h = plot(A{i}(:,1), A{i}(:,4), lw, 1.5);
  end
  K = A{i}(1,2);
  labels{i} = ['k = ' num2str(K)];
  colours{i} = get(h, 'color');
end
set(gca, 'fontsize', 14)

for i=1:5
  h = plot(A{i}(:,1), A{i}(:,3), 'k--', lw, 1);
  K = A{i}(1,2);
  set(h, 'color', colours{i})
  labels{end+1} = ['k = ' num2str(K) ' quad'];
end

N = A{i}(:,1);
xlabel('N')
ylabel('error')
if d == 2
  plot(N, 0.004*N.^(-0.5), 'k:', lw, 1.5)
  plot(N, 0.008*(log(N)./N), 'k--', lw, 1.5)
  plot(N, 0.033*N.^(-1), 'k:', lw, 1.5)
  labels{end+1} = 'N^{-1/2}';
  labels{end+1} = 'log(N)/N';
  labels{end+1} = 'N^{-1}';
  legend(labels, 'location', 'eastoutside')
  ylim([1e-7 4e-3])
  xlim([1e2 4e6])
  title('d = 2')
elseif d == 3
  plot(N, 0.07*N.^(-0.5), 'k:', lw, 1.5)
  h = plot(N, 0.05*(log(N)./N).^(2/3), 'k--', lw, 1.5)
  plot(N, 0.66*N.^(-1), 'k:', lw, 1.5)
  %labels{end+1} = 'N^{-1/2}';
  %labels{end+1} = '(log(N)/N)^{2/3}';
  %labels{end+1} = 'N^{-1}';
  %legend(labels, 'location', 'eastoutside')
  legend(h, {'(log(N)/N)^{2/3}'});
  ylim([1e-6 2e-2])
  xlim([1e2 4e6])
  title('d = 3')
elseif d == 4
  plot(N, 0.06*N.^(-0.5), 'k:', lw, 1.5)
  h = plot(N, 0.03*(log(N)./N).^(2/4), 'k--', lw, 1.5)
  plot(N, 0.6*N.^(-1), 'k:', lw, 1.5)
  %labels{end+1} = 'N^{-1/2}';
  %labels{end+1} = '(log(N)/N)^{1/2}';
  %labels{end+1} = 'N^{-1}';
  legend(h, {'(log(N)/N)^{1/2}'});
  ylim([1e-5 2e-2])
  xlim([1e2 4e6])
  title('d = 4')
else
  error('invalid dimension')
end
%xlim([100 300000])
b = ['line_integral_uniform_d' num2str(d)']
print(b, '-dpng')
print(b, '-depsc2')
assert(~ system(['epstopdf ' b '.eps']))
end
