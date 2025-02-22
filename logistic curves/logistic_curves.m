resolution = [510, 400];
figure
set(gcf, Position=[0, 0, resolution], Renderer='painter')

n_i = 400; % number of frames
filename = 'logistic curves dark pulse hue annotated 510x400';
savegif = true;
savevid = true;
if savegif || savevid
	F(n_i) = struct('cdata',[],'colormap',[]);
	if savevid
		V = VideoWriter([filename, '.mp4'], 'MPEG-4');
		set(V, FrameRate=24, Quality=100);
		open(V)
	end
end

% watermark
annotation('TextBox', ...
	String='@jaketrobert1000', ...
	Position=[1, 0, 0, 0], ...
	HorizontalAlignment='right', ...
	VerticalAlignment='bottom', ...
	FitBoxToText='on', ...
	BackgroundColor='k', ...
	Color='w', ...
	LineStyle='none',...
	Margin=3)

% caption
T_caption = annotation('TextBox', ...
	Interpreter='latex', ...
	FontSize=16, ...
	Position=[0 + 7/resolution(1), 1 - 7/resolution(2), 0, 0], ...
	HorizontalAlignment='left', ...
	VerticalAlignment='top', ...
	FitBoxToText='on', ...
	BackgroundColor='k', ...
	Color='w', ...
	EdgeColor='w',...
	Margin=5);

load c_logistic_curves.mat

for i=1:n_i
	L = 15; % carrying capacity
	k = 0.3; % growth coefficient
	n_curves = 144;
	N0 = linspace(-40, 0, n_curves) + i*55/n_i; % initial population
	% N0 = linspace(-40, 0, n_curves) + 20;
	
	t = linspace(-20, 30, 200)';
	N = L*N0.*exp(1).^(k*t)./(L + N0.*(exp(1).^(k*t) - 1));
	N(find(abs(N) > 10^2)) = NaN;

	c = interp1(linspace(0, 1, size(c_logistic_curves, 1)), ...
		c_logistic_curves, ...
		erf(10*i/n_i)*sin(2*linspace(0, 1, n_curves)*pi - 10*i/n_i*pi).^2);
	c_hsv = rgb2hsv(c);
	c_hsv(:, 1) = clip(c_hsv(:, 1) + linspace(0, 1, n_curves)', 0, 1);
	c = hsv2rgb(c_hsv);
	% c = turbo(n_curves);
	
	P = plot(t, N, LineWidth=1);
	xlim([-20, 30]), ylim([-15, 30])
	set(gca, Position=[0, 0, 1, 1], ...
		xAxisLocation='origin', ...
		yAxisLocation='origin', ...
		xTick=-15:5:25, ...
		yTick=-10:5:25, ...
		xColor='w', ...
		yColor='w', ...
		Layer='top', ...
		Color='k')
	colororder(c)
	box off

	T_caption.String = ['$N(t)=\frac{15N_0e^{0.3t}}{15 + N_0(e^{0.3t}-1)} \;\; \textbf{for} \; ', num2str(N0(1), '%.1f'),' < N_0 < ', num2str(N0(n_curves), '%.1f'),'$']

	drawnow

	if savegif || savevid
		F(i) = getframe(gcf);
		if savegif
			[ind, map] = rgb2ind(frame2im(F(i)), 32);
			if i == 1
				imwrite(ind, map, [filename, '.gif'], 'gif', DelayTime=1/24, LoopCount=Inf)
			else
				imwrite(ind, map, [filename, '.gif'], 'gif', DelayTime=1/24, WriteMode='append')
			end
		end
		if savevid
			writeVideo(V, F(i))
		end
	end
end

if savegif || savevid
	if savevid
		close(V)
	end
	fig = figure;
	set(gcf, Position=[0, 0, resolution])
	movie(fig, F, 1)
end