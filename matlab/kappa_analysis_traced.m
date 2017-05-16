%function kappa_analysis_traced(input_file)
% function for analysing CSV files of kappa values and to display tracing result
% of the signatures in the CSS

% 0 = Position Color (Rainbow)
% 1 = Random Color
% 2 = Red/Green
colorOption = 0;

% define some colors
sharkblue = [0.37,0.55,0.83];

% read file
data_in = csvread('/home/max/gits/css-new/build/Data/kappaAppr.csv');
%data_in = csvread(input_file);

% first column contains sigma values
sigma_values = data_in(1:end, 1);

% number of columns in data_in
num_cols_in = size(data_in, 2);

% index vector (num_cols_in-2 because first value is sigma value)
u_index = 0:1:(num_cols_in-2);

% new figure
fig2 = figure;
grid on;
hold on;
pos_fig2 = [0 0 1200 900]; % sets the size of the figure in px
set(fig2, 'Position', pos_fig2);
view(3); % view(az,el) - view(3) means az = -37.5, el = 30
% set axis labels
xlabel('u (arc length)');
ylabel('\sigma (scale)');
zlabel('\kappa (curvature)');

% plot every selected line of kappa as one graph
sigmaStart = 1; % number of start value (default 40)
sigmaStep = 1; % default 10 (only every 10th value)
for ind = sigmaStart:sigmaStep:length(sigma_values)
	current_sigma = sigma_values(ind) * ones(1, num_cols_in-1);
	current_kappa = data_in(ind, 2:end); % row ind, col 2 to end
	plot3(u_index, current_sigma, current_kappa, 'Color', sharkblue);
end

% read minima
minima_in = csvread('/home/max/gits/css-new/build/Data/tracedMin.csv');

% read how many values the traced signatures have
minima_nu_values = minima_in(:,1);

% save color value for every traced signature
minima_rgb_values = minima_in(:,2:4);

% the 4th value in a row is the first u-value
% ind1 goes through the lines
for ind1 = 1:1:length(minima_nu_values)
	% break if no entries
	if minima_nu_values(1,1) == 0
		break;
	end

	% color value
	if colorOption == 0
		rgb_val = [minima_rgb_values(ind1, 1)/255, minima_rgb_values(ind1, 2)/255, minima_rgb_values(ind1, 3)/255];
	elseif colorOption == 1
		rgb_val = rand(3,1);
	else
		rgb_val = [1, 0, 0];
	end

	% ind2 goes through columns
	for ind2 = 5:3:(3*minima_nu_values(ind1)+2)
		plot3(minima_in(ind1, ind2), minima_in(ind1, ind2+1), minima_in(ind1, ind2+2), 'v', 'Color', rgb_val);
	end
	
	% print number of min at bottom or side
	strNumber = num2str(ind1-1);
	os = 3*minima_nu_values(ind1)+2;
	%text(minima_in(ind1, 5)-7, minima_in(ind1, 6)-4, minima_in(ind1, 7), strNumber, 'Color', rgb_val);
	text(minima_in(ind1, os)+10, minima_in(ind1, os+1), minima_in(ind1, os+2), strNumber, 'Color', rgb_val);
end

% read maxima
maxima_in = csvread('/home/max/gits/css-new/build/Data/tracedMax.csv');

% read how many values the traced signatures have
maxima_nu_values = maxima_in(:,1);

% save color value for every traced signature
maxima_rgb_values = maxima_in(:,2:4);

% the 4th value in a row is the first u-value
for ind1 = 1:1:length(maxima_nu_values)
	% break if no entries
	if maxima_nu_values(1,1) == 0
		break;
	end

	% color value
	if colorOption == 0
		rgb_val = [maxima_rgb_values(ind1, 1)/255, maxima_rgb_values(ind1, 2)/255, maxima_rgb_values(ind1, 3)/255];
	elseif colorOption == 1
		rgb_val = rand(3,1);
	else
		rgb_val = [0, 1, 0];
	end
	
	for ind2 = 5:3:(3*maxima_nu_values(ind1)+2)
		plot3(maxima_in(ind1, ind2), maxima_in(ind1, ind2+1), maxima_in(ind1, ind2+2), '^', 'Color', rgb_val);
	end
	
	% print number of max at bottom or side
	strNumber = num2str(ind1-1);
	os = 3*maxima_nu_values(ind1)+2;
	%text(maxima_in(ind1, 5)-7, maxima_in(ind1, 6)-4, maxima_in(ind1, 7), strNumber, 'Color', rgb_val);
	text(maxima_in(ind1, os)+10, maxima_in(ind1, os+1), maxima_in(ind1, os+2), strNumber, 'Color', rgb_val);
end

% ('MarkerFaceColor', rgb_val) for filling the markers



