% init
figure
resolution = [360, 360];
margin_x = 0.05;
margin_y = 0.5*(1 - resolution(1)/resolution(2)*(1 - 2*margin_x));
L_x = 1 - 2*margin_x;
L_y = 1 - 2*margin_y;
set(gcf, Position=[0, 0, resolution], Renderer='painter', Color='k')
set(gca, Position=[margin_x, margin_y, L_x, L_y])

T = annotation('TextBox', String='@jaketrobert1000', Position=[0, 0, 0, 0], VerticalAlignment='bottom', FitBoxToText='on', LineStyle='none', BackgroundColor='k', Color='w')

load c.mat

% main
lim = 320;
% n0 = 12;
n0 = 8; % # of periods
n1 = (1:12); % # of periods per period
n_frames = 24*8;
n_vertices = n0*10^2;

savevideo = true;
filename = 'swirls inverted 360x360';
if savevideo
	F(n_frames) = struct('cdata', [], 'colormap', []);
	V = VideoWriter([filename, '.mp4'], 'MPEG-4');
	set(V, FrameRate=24, Quality=100);
	open(V)
end

for i = 1:n_frames
	a = power(i/n_frames, 1/2)*(1 - 0.4*sin(4*pi*i/n_frames)^2);
	c = [interp1(linspace(0, 1, size(c_red, 1)), flip(c_cyan, 1), a*linspace(0, 1, size(n1, 2)).^4); ...
		interp1(linspace(0, 1, size(c_cyan, 1)), flip(c_red, 1), a*[linspace(0, 1, size(n1, 2)).^4, 1])];

	t = linspace(0, 2*pi, n_vertices+1)';

	scale0 = interp1([0, 1], [-100, -32], power(i/n_frames, 1/4));
	scale1 = interp1([0, 1], [60, 212-60], i/n_frames) * (1 + 0.05*sin(2*t + 6*2*pi*i/n_frames)); % 0 --> 225 (default 212)
	% scale1 = interp1([0, 1], [40, 212-60], i/n_frames) * (1 + 0.05*sin(2*t + 6*2*pi*i/n_frames)); % 0 --> 225 (default 212)
	
	r1 = interp1([0, 1], [0, 60], i/n_frames); % inner radius
	r2 = interp1([0, 1], [310, 300], i/n_frames); % outer radius
	% phi = 0.01*sin(linspace(0, pi, size(n1, 2)) + 0.5*2*pi*i/n_frames);
	% phi = pi/n0;
	phi = 0;
	
	weight = (1 - abs(cos(n0*t/2))).^2;
	% r = weight.*(scale1 - scale0*abs(cos(n0*n1.*t)));
	r = weight.*(scale1 - scale0*cos(n0*n1.*t));
	theta = t - weight*scale0.*sin(n0*n1.*t)./(0.2*r*n0);
	
	x1 = (r1 + r).*cos(theta + pi/n0 + phi);
	y1 = (r1 + r).*sin(theta + pi/n0 + phi);
	x2 = (r2 - r).*cos(theta + phi);
	y2 = (r2 - r).*sin(theta + phi);
	x3 = 1.01*r2*cos(linspace(0, 2*pi, 4*n0+1));
	y3 = 1.01*r2*sin(linspace(0, 2*pi, 4*n0+1));
	
	P = plot(x1, y1, x2, y2, x3, y3);
	colororder(c)
	xlim([-lim, lim]), ylim([-lim, lim])
	set(P, LineWidth=1, LineJoin='chamfer');
	% P(size(n1, 2)).LineWidth = 1;
	% P(2*size(n1, 2)).LineWidth = 1;

	set(gca, xTick=[], yTick=[], Clipping='off')
	set(gca, Color='k', xColor='k', yColor='k')

	drawnow
	if savevideo
		F(i) = getframe(gcf);
		[ind, map] = rgb2ind(frame2im(F(i)), 16);
		if i == 1
			imwrite(ind, map, [filename, '.gif'], 'gif', DelayTime=1/24, LoopCount=Inf)
		else
			imwrite(ind, map, [filename, '.gif'], 'gif', DelayTime=1/24, WriteMode='append')
		end
		writeVideo(V, F(i))
	end
end

if savevideo
	close(V)

	% figure
	% set(gcf, Position=[0, 0, resolution])
	% movie(F, 1)
end