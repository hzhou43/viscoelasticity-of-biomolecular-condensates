function output = dropletrupture(Filename,datacolor1,n,tstart,tend,t0,orientation,radius)

% Input the filename (as in 'rupturefile.h5', the color for plotting,
% the moving average window size n (default 5000),
% the limits of the rupture period from time = tstart to tend, and
% the point t0 in the x axis after which the force begins to show change.
% Orientation is 'lt' or 'rt'. A downward dipping curve for the recorded force with trap1 will
% always have a left orientation and a upward rising curve will have a
% right orientation. radius is the bead radius.

F1x = h5read(Filename,'/Force HF/Force 1x');

timeindex =[0:1:length(F1x)-1];
timeindex = timeindex';
t = timeindex./78125;

%moving average
F1x = movmean(F1x,n);

figure
hold
z = F1x(t>t0);
m2 = F1x(1);

if (orientation=='rt') 
    yline(max(z)-m2,'b');
    m1 = max(z);
elseif (orientation=='lt')
    yline(min(z)-m2,'b');
    m1 = min(z);
end
yline(F1x(1)-m2,'r');

xlim([0 tend-tstart])
plot(t,F1x-m2,datacolor1)
xlabel('Time (s)')
grid on

ST = abs((m1-m2)/(2*pi*radius*1.1));

output = ['The surface tension is ',num2str(ST),' pN/um'];

end
