%% init
resolution = [
	720, 480;
];
margin_x = 0.06;
% margin_y = 0.5*(1 - resolution(1)/resolution(2)*(1 - 2*margin_x));
margin_y = margin_x*resolution(1)/resolution(2);
L_x = 1 - 2*margin_x;
L_y = 1 - 2*margin_y;

figure
set(gcf, ...
	Position = [0, 0, resolution(1, :)], ...
	Color = 'w', ...
	Renderer = 'painter' ...
)
set(gca, ...
	Position = [margin_x, margin_y, L_x, L_y] ...
)

% init video
filename = ['rc voltage alt ', num2str(resolution(1)), 'x', num2str(resolution(2))];
savegif = true;
savevid = true;
if savevid
	vid = VideoWriter([filename, '.mp4'], 'MPEG-4');
	set(vid, FrameRate=24, Quality=50)
	open(vid)
end

% watermark
annotation('TextBox', ...
	Position = [1, 0, 0, 0], ...
	String = '@jaketrobert1000', ...
	FitBoxToText = 'on', ...
	HorizontalAlignment = 'right', ...
	VerticalAlignment = 'bottom', ...
	BackgroundColor = get(gcf, 'color'), ...
	Color = 'k', ...
	LineStyle = 'none' ...
)

%% main
a = 4.0; % linear frequency constant
b = 0.5; % exponential frequency constant

% extended: 24000, normal: 500
x0 = lambertw(0, 1.0*pi*(0:24000)*b/a)/b; % sample points
n_x0 = size(x0, 2);
V_pp = @(x) 0.25*exp(1).^(0.60*x); % peak-to-peak voltage function
V0 = 0.5*V_pp(x0).*cos(a*x0.*exp(1).^(b*x0)); % square wave voltage signal

v = 1/24; % sweep velocity (sec/frame)
lim_x = 4; % x-limits = [0, lim_x] + offset
n_vertices = 16; % n of vertices per voltage point
n_i = floor((x0(end) - lim_x)/v); fprintf('n of frames: %d \n', n_i) % n of frames
P_V = stairs(x0, V0, LineWidth=1, Color='k'); % voltage plot
line_x = xline(0, LineStyle=':', LineWidth=1, Layer='bottom', Color='k'); % sweep line
% line_y = yline(0, Layer='top', Color=); % current voltage line

T = logspace(3.1053, -0.5000, 20); % time constant
% T = 10^3.1053;
n_T = size(T, 2);
V_prev = zeros(1, n_T);
k_prev = 0;

if savegif || savevid
	frame(n_i) = struct('cdata', [], 'colormap', []);
end

% annotation
T1 = annotation('TextBox', Position=[margin_x, margin_y/2 - 3/resolution(2), 0, 0], LineStyle='none');
T2 = annotation('TextBox', Position=[margin_x+0.18*1, margin_y/2, 0, 0],  LineStyle='none');
T3 = annotation('TextBox', Position=[0.5 + 6/resolution(1), 0.5 - 16/resolution(2), 0, 0], FaceAlpha=0.85, LineStyle='none');
set([T1, T2, T3], ...
	FitBoxToText = 'on', ...
	HorizontalAlignment = 'left', ...
	VerticalAlignment = 'middle', ...
	BackgroundColor = get(gcf, 'color'), ...
	Color = 'k', ...
	FontWeight = 'bold', ...
	Margin = 2 ...
)

title('Voltage vs. time (in the evil RC circuit)')

% loop
x = zeros(n_x0*n_vertices, 1)*NaN;
V = zeros(n_x0*n_vertices, n_T)*NaN;
hold on
P_lgd = plot(x, V, LineWidth=3);
for i = 1:n_i
	x_i = (i - 1)*v;
	k = find(x0 > x_i, 1); % index of last voltage change
	if k > 2 && k ~= k_prev
		if k > 3
			delete(P_prev)
		end
		for j = k_prev:(k - 1)
			x_j = (linspace(x0(j - 1), x0(j), n_vertices)');
			x((1:n_vertices) + (j - 1)*n_vertices) = x_j;
			V((1:n_vertices) + (j - 1)*n_vertices, :) = (V_prev - V0(j - 1)).*exp(1).^(T.*(x0(j - 1) - x_j)) + V0(j - 1);
			V_prev = V(n_vertices*(j), :);
		end
		P_prev = plot(x, V, LineWidth=1);
	end

	x0_k = x0(k - 1);
	V0_k = V0(k - 1);
	k_range = (1:n_vertices) + (k - 2)*n_vertices;
	x_k = (linspace(x0_k, x_i, n_vertices)');
	V_k = (V_prev - V0_k).*exp(1).^(T.*(x0_k - x_k)) + V0_k;
	
	
	P = plot(x_k, V_k, 'k', LineWidth=1);
	P_marker = plot(x_i, V_k(end, :), '.', MarkerSize=8);
	P_V_line = plot([x0_k, x_i], [V0_k, V0_k], 'k-', LineWidth=2);
	P_V_marker = plot(x_i, V0_k, 'k*', MarkerSize=6);
	box on
	line_x.Value = x_i; % update sweep line
	ax = set(gca, ...
		xAxisLocation = 'origin', ...
		yAxisLocation = 'origin', ...
		xTick = 0, ...
		xTickLabel = '0       ', ...
		yTick = [], ...
		Layer = 'bottom', ...
		Color = 'w', ...
		xColor = 'k', ...
		yColor = 'k' ...
	);
	% xlabel('Time (s)');
	% ylabel('Voltage (V)');
	xlim([0, lim_x] + x_i - 0.5*lim_x)
	ylim([-1, 1]*0.40*V_pp(0.5*lim_x + x_i)/2)

	% color
	c0 = hex2rgb('#99b31b');
	c_phase = 0; % 0.001*i;
	for j = 1:n_T
		c_j = hsv2rgb(mod(rgb2hsv(c0) + [c_phase + 1.00*j/n_T, 0, 0], 1.0));
		P(j).Color = c_j;
		P_marker(j).Color = c_j;
		if k > 2
			P_prev(j).Color = c_j;
		end
		if i == 1
			P_lgd(j).Color = c_j;
		end
	end

	% update annotation
	% f = a*exp(1).^(b*x_i)/(2*pi);
	f = a*exp(1).^(b*x_i)*(1 + b*x_i)/(2*pi);
	str_f = sprintf('f = %.3f Hz', f);
	str_x = sprintf('t = %.3f s', x_i);
	str_V = sprintf('V_{0} = %.3f V_{pp}', 2*abs(V0_k));
	if i == 1
		ind_lgd = 1:3:size(P_lgd, 1);
		str_lgd = cell(size(ind_lgd, 2), 1);
		for j = 1:size(ind_lgd, 2)
			str_lgd{j} = sprintf('T = %.1e Ω^{-1}F^{-1}', T(ind_lgd(j)));
		end
		lgd = legend(P_lgd(ind_lgd), str_lgd, AutoUpdate='off', IconColumnWidth=10);
		lgd.Title.String = 'Capacitor time constant';
		lgd.Direction = 'reverse';
	end
	T1.String = str_V;
	T2.String = str_f;
	T3.String = str_x;
	
	drawnow

	% update
	k_prev = k;
	if savegif || savevid
		frame(i) = getframe(gcf);
		[ind, map] = rgb2ind(frame2im(frame(i)), 32);
		if savegif
			if i == 1
				imwrite(ind, map, [filename, '.gif'], 'gif', DelayTime=1/24, LoopCount=Inf)
			else
				imwrite(ind, map, [filename, '.gif'], 'gif', DelayTime=1/24, WriteMode='append')
			end
		end
		if savevid
			writeVideo(vid, frame(i));
		end
	end
	
	% clear
	if i ~= n_i
		delete(P)
		delete(P_marker)
		delete(P_V_line)
		delete(P_V_marker)
	end
end
%% end
if savevid
	close(vid)
end