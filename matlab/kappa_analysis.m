%function kappa_analysis(input_file)
% function for analysing CSV files of kappa values

% define some colors
sharkblue = [0.37,0.55,0.83];

% read file
%data_in = csvread('kappa.csv');
data_in = csvread('/home/max/gits/css-new/build/Data/kappaAppr.csv');
%data_in = csvread(input_file);

% first column contains sigma values
sigma_values = data_in(1:end, 1);

% number of columns in data_in
num_cols_in = size(data_in, 2);

% index vector
u_index = 0:1:(num_cols_in-2);

% new figure
fig1 = figure;
grid on;
hold on;
pos_fig1 = [0 0 1200 900]; % sets the size of the figure in px
set(fig1, 'Position', pos_fig1);
view(3); % view(az,el) - view(3) means az = –37.5, el = 30
% set axis labels
xlabel('u (arc length)');
ylabel('\sigma (scale)');
zlabel('\kappa (curvature)');

% plot every selected line of kappa as one graph
sigmaStart = 1; % number of start value (default 40)
sigmaStep = 1; % default 10
for ind = sigmaStart:sigmaStep:length(sigma_values)
	current_sigma = sigma_values(ind) * ones(1, num_cols_in-1);
	current_kappa = data_in(ind, 2:end); % row ind, col 2 to end
	plot3(u_index, current_sigma, current_kappa, 'Color', sharkblue);
end

% read minima values
minima_in = csvread('/home/max/gits/css-new/build/Data/min.csv');

% plot every selected line of minima as one graph
% note that matlab fills in zeros
for ind = sigmaStart:sigmaStep:length(sigma_values)
	current_number_min = minima_in (ind, 2);
	current_min_line = minima_in (ind, 3:(2+current_number_min));

	% find values which are not zero
	%current_min_line = minima_in (ind, 2:end);
	%current_min_line(current_min_line == 0) = []; % = [] means to delete the values if true
	
	% get corresponding kappa values
	current_kappa = data_in(ind, 2:end); % row ind, col 2 to end
	current_kappa = current_kappa(current_min_line+1); % matlab index starts at 1
	
	% get current sigma value and create corresponding vector
	current_sigma = sigma_values(ind) * ones(1, length(current_min_line));
	
	% plot
	plot3(current_min_line, current_sigma, current_kappa,'r*');
end

% read maxima values
maxima_in = csvread('/home/max/gits/css-new/build/Data/max.csv');

for ind = sigmaStart:sigmaStep:length(sigma_values)
	current_number_max = maxima_in (ind, 2);
	current_max_line = maxima_in (ind, 3:(2+current_number_max));

	% find values which are not zero
	%current_max_line = maxima_in (ind, 2:end);
	%current_max_line(current_max_line == 0) = []; % = [] means to delete the values if true
	
	% get corresponding kappa values
	current_kappa = data_in(ind, 2:end); % row ind, col 2 to end
	current_kappa = current_kappa(current_max_line+1); % matlab index starts at 1
	
	% get current sigma value and create corresponding vector
	current_sigma = sigma_values(ind) * ones(1, length(current_max_line));

	% plot
	plot3(current_max_line, current_sigma, current_kappa,'g*');
end

% ************ Visualization ************
% rotate view to 2D
rotatePlot = 0;
saveImageSequence = 0;

if rotatePlot == 1
	pause(10); % show normal plot for n seconds

	az = -37.5;
	el = 30;

	% save first frame
	if saveImageSequence == 1
		% save as png
		print('seq/00','-dpng')
	end

	for ind = 1:1:60
		az = az+0.625;
		el = el+1;
		
		pause(0.1);
		view(az, el);
		
		if saveImageSequence == 1
			fileName = num2str(ind);
			fileName = strcat('seq/',fileName);
			print(fileName,'-dpng');
		end
	end

	pause(3);
	% rotate back again
	for ind = 1:1:60
		az = az-0.625;
		el = el-1;
		
		pause(0.1);
		view(az, el);
	end
end
