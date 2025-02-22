figure
resolution = [480, 480];
margin_x = 0.1;
margin_y = 0.5*(1 - resolution(1)/resolution(2)*(1 - 2*margin_x));
L_x = 1 - 2*margin_x;
L_y = 1 - 2*margin_y;
set(gcf, Position=[0, 0, resolution], ...
	Color='#000', ...
	Renderer='painter')
set(gca, Position=[margin_x, margin_y, L_x, L_y])

filename = ['love equation dark ', num2str(resolution(1)), 'x', num2str(resolution(2))];
disp(filename)
savegif = true;
savevid = true;
ni = 200;
if savegif || savevid
	F(ni) = struct('cdata',[],'colormap',[]);
	if savevid
		V = VideoWriter([filename, '.mp4'], 'MPEG-4');
		set(V, FrameRate=24, Quality=100);
		open(V)
	end
end

annotation('TextBox', ...
	String='@jaketrobert1000', ...
	Position=[1, 0, 0, 0], ...
	Color='#fff', ...
	FitBoxToText='on', ...
	HorizontalAlignment='right', ...
	VerticalAlignment='bottom', ...
	BackgroundColor=get(gcf, 'color'), ...
	Margin=2, ...
	LineStyle='none')

% latex

% annotation('TextBox', ...
% 	String=[
% 		"$r(t)= 6.4 + 3.0cos(\omega t) + 2.4(1 - |cos(0.5\omega t)|)^2$", ...
% 		"$\theta(t)= t + 2.6sin(\omega t) \cdot r(t)^{-1}$"
% 	], ...
% 	Interpreter='latex', ...
% 	FontSize=12, ...
% 	Position=[0.5, 1 - 0.67*margin_y, 0, 0], ...
% 	FitBoxToText='on', ...
% 	HorizontalAlignment='center', ...
% 	VerticalAlignment='middle', ...
% 	BackgroundColor=get(gcf, 'color'), ...
% 	Margin=2, ...
% 	LineStyle='none')


% tex

annotation('TextBox', ...
	String=[
		"r(t)= 6.4 + 3.0cos(ωt) + 2.4(1 - |cos(0.5ωt)|)^2", ...
		"θ(t)= t + 2.6sin(t)*r(t)^{-1}"
	], ...
	Interpreter='tex', ...
	Position=[0.5, 1 - 0.5*margin_y, 0, 0], ...
	FitBoxToText='on', ...
	HorizontalAlignment='center', ...
	VerticalAlignment='middle', ...
	BackgroundColor=get(gcf, 'color'), ...
	Margin=2, ...
	LineStyle='none')


% T = annotation('TextBox', ...
% 	String=[
% 		"$u(t)=6.75A \Bigg(( 1 - \sqrt{1 - mod\Big(\frac{(-t+\phi)}{\pi}, 1\Big)^2} \; \Bigg)^2 \cdot \sqrt{1 - mod(\frac{(-t+\phi)}{\pi}, 1)} \;\; \textbf{for:}$", ...
% 		"$\phi=1, \;\; -1.0 < A < 8.0$"
% 	], ...
% 	Interpreter='latex', ...
% 	FontSize=10, ...
% 	Position=[0.5, 0.5*margin_y, 0, 0], ...
% 	FitBoxToText='on', ...
% 	HorizontalAlignment='center', ...
% 	VerticalAlignment='middle', ...
% 	BackgroundColor=get(gcf, 'color'), ...
% 	Margin=2, ...
% 	LineStyle='none');

% load c0.mat

w1 = 8.0; % frequency
w2 = -4.0; % pulse frequency
lim = 10;
for i = 1:ni
	phi = 2*pi*i/ni;

	t = linspace(0, 2*pi, w1*10^2 + 1)';
	% pulse = 1.0*cos(t/2 + phi).^2; % cos^2 function
	% pulse = 1.0*mod((t - phi)/(2*pi), 1).^3; % sawtooth function
	pulse = (1 - sqrt(1 - mod(-(t + phi)/(1.0*pi), 1).^2)).^2.*sqrt(1 - mod(-(t + phi)/(1.0*pi), 1))/0.14815; % pulse function
	pulse_phi = -0.8*2*pi*i/100;
	r1 = 6.4 + 0.1*sin(phi) + 3.0*cos(w1*t) + 2.40*(1 - abs(cos(w1*t/2))).^2 + pulse.*(1.0).*linspace(12.0, -1.0, 24);
	r2 = 3.5 - 0.1*sin(phi) - 1.67*cos(w1*(t + pi/w1)) - 2.10*(1 - abs(cos(w1*(t + pi/w1)/2))).^2 - 0.240*pulse.*(1.0).*linspace(10.0, -1.0, 12);
	% theta = t + 2.60*sin(w1*t)./r + pulse.*sin(w2*(t + pulse_phi))./r;
	theta1 = t + 2.60*sin(w1*t)./r1;
	theta2 = t + 1.60*sin(w1*(t + pi/w1))./r2;
	
	x1 = r1.*cos(theta1);
	y1 = r1.*sin(theta1);
	x2 = r2.*cos(theta2);
	y2 = r2.*sin(theta2);
	
	xc1 = sin(1*phi)^2;
	xc2 = sin(1*phi - 0.25*pi)^2;

	% c1_r = interp2(0:1, linspace(0, 1, 256), [c0_1(:,1), c0_2(:,1)], xc1, linspace(0, 1, size(r1, 2)));
	% c1_g = interp2(0:1, linspace(0, 1, 256), [c0_1(:,2), c0_2(:,2)], xc1, linspace(0, 1, size(r1, 2))); 
	% c1_b = interp2(0:1, linspace(0, 1, 256), [c0_1(:,3), c0_2(:,3)], xc1, linspace(0, 1, size(r1, 2)));
	c1_r = interp2(0:1, linspace(0, 1, 256), [c0_1d(:,1), c0_2d(:,1)], xc1, linspace(0, 1, size(r1, 2)));
	c1_g = interp2(0:1, linspace(0, 1, 256), [c0_1d(:,2), c0_2d(:,2)], xc1, linspace(0, 1, size(r1, 2))); 
	c1_b = interp2(0:1, linspace(0, 1, 256), [c0_1d(:,3), c0_2d(:,3)], xc1, linspace(0, 1, size(r1, 2)));
	c1 = [c1_r, c1_g, c1_b];
	
	% c2_r = interp2(0:1, linspace(0, 1, 256), [c0_1(:,1), c0_2(:,1)], xc2, linspace(0, 1, size(r2, 2)));
	% c2_g = interp2(0:1, linspace(0, 1, 256), [c0_1(:,2), c0_2(:,2)], xc2, linspace(0, 1, size(r2, 2))); 
	% c2_b = interp2(0:1, linspace(0, 1, 256), [c0_1(:,3), c0_2(:,3)], xc2, linspace(0, 1, size(r2, 2)));
	c2_r = interp2(0:1, linspace(0, 1, 256), [c0_1d(:,1), c0_2d(:,1)], xc2, linspace(0, 1, size(r2, 2)));
	c2_g = interp2(0:1, linspace(0, 1, 256), [c0_1d(:,2), c0_2d(:,2)], xc2, linspace(0, 1, size(r2, 2))); 
	c2_b = interp2(0:1, linspace(0, 1, 256), [c0_1d(:,3), c0_2d(:,3)], xc2, linspace(0, 1, size(r2, 2)));
	c2 = [c2_r, c2_g, c2_b];
	
	P1 = plot(x1, y1, LineWidth=1); hold on
	P2 = plot(x2, y2, LineWidth=1); hold off
	P1(end).LineWidth = 1;
	P2(end).LineWidth = 1;
	xlim([-lim, lim]), ylim([-lim, lim])
	% set(gca, xTick=-10:5:10, yTick=-10:5:10)
	set(gca, xTick=[], yTick=[])
	set(gca, ...
		Color='#000', ...
		xColor='#fff', ...
		yColor='#fff', ...
		Layer='top')
	colororder([c1; c2])

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