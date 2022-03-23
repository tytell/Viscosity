function scallopVideo(outname,inertia)

fps = 15;
freq = 0.8;
dt = 1/fps;
t = 0:dt:4/freq;
ang = (50*cos(2*pi*freq*t) + 70)*pi/180;     % 120 down to 20 deg
angvel = -50*2*pi*sin(2*pi*freq*t)*pi/180;

linpos = 0.5*(-100*pi/180/(2*pi*freq)*cos(2*pi*freq*t) + 0.5*inertia*t.^2);

len = 1;
thick = 0.1;

len2 = len/2;
thick2 = thick/2;
ang2 = ang/2;

x0 = [len2 -len2 -len2 len2];
y0 = [-thick2 -thick2 thick2 thick2];

fig = figure(1);
clf;
set(fig,'Units','pixels','WindowStyle','normal','Color','w');
pos = get(fig,'Position');
pos(3) = 640;
pos(4) = 480;
set(fig,'Position',pos);

ax = axes('Parent',fig,'Position',[0 0 1 1]);
hold on;

x1 = x0*cos(ang2(1)) + y0*sin(ang2(1));
y1 = -x0*sin(ang2(1)) + y0*cos(ang2(1));
y1 = y1 - y1(1);
x1 = x1 - x1(1);

xx = x1([1 2 3 4 4 3 2]);
yy = cat(2,y1([1 2 3 4]),-y1([4 3 2]));

hs = fill(xx,yy,[0.7 0.7 0.7],'EdgeColor','none');
axis equal;
yl = 3*480/640;
axis(ax, [-1 2 -yl/2 yl/2]);
hold off;
axis off;

if (~isempty(outname))
    vid = VideoWriter(outname);
    vid.FrameRate = fps;
    open(vid);
end

for i = 1:length(t)
    x1 = x0*cos(ang2(i)) + y0*sin(ang2(i));
    y1 = -x0*sin(ang2(i)) + y0*cos(ang2(i));
    y1 = y1 - y1(1);
    x1 = x1 - x1(1);

    xx = x1([1 2 3 4 4 3 2]);
    yy = cat(2,y1([1 2 3 4]),-y1([4 3 2]));
    
    set(hs,'XData',xx+linpos(i), 'YData',yy);
    drawnow;
    
    if (~isempty(outname))
        frame = getframe(ax);
        writeVideo(vid,frame);
    else
        drawnow;
        pause;
    end
end    
if ~isempty(outname)
    close(vid);
end




    