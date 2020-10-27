%% Estimate densities for all years and countries
addpath("~/your/working/directory/Others/hades/codes")

load('~/your/working/directory/Data/HMD/age_counts.mat') % contains variables ctry, x, tms and yrs

mkdir('~/your/working/directory/Data/HMD/dens_est')

ind = find(x > 20);  % restrict Mortality to ages greater than 20
x = x(ind);
for i = 1:length(ctry)
	disp(i)
	yr = eval(strcat('yrs.', ctry{i}));
	tm = eval(strcat('tms.', ctry{i})); tm = tm(ind, :);
	dns = zeros(length(x), length(yr));

	for j = 1:length(yr)

		y = cell2mat(arrayfun(@(k) x(k)*ones(1, tm(k, j)), 1:length(x), 'UniformOutput', false));
		tmp = hist2fx(y, 1, 2, 'gauss', x);
		tmp = max(tmp.dens, 0);
		dns(:,j) = tmp/trapz(x, tmp);

	end

	save(strcat('~/your/working/directory/Data/HMD/dens_est/', ctry{i}, '_est'), 'dns', 'x', 'yr')

end
