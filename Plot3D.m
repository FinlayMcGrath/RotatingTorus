%set initial variables
h = 0.05;
t = 16;
nmax = t / h;
a = 1;
c = 4;
M = 3.14;
%
wx = zeros(1,nmax);
wy = zeros(1,nmax);
wz = zeros(1,nmax);

wdotx = zeros(1,nmax);
wdoty = zeros(1,nmax);
wdotz = zeros(1,nmax);

%impulses
I1 = M * (((5/8) * a^2) + ((1 / 2) * c^2));
I2 = M * (((5/8) * a^2) + ((1 / 2) * c^2));
I3 = M * (((3/4) * a^2) + c^2);
I = [I1, 0, 0, 
    0, I2, 0, 
    0, 0, I3];

y1 = (I3 - I2) / I1;
y2 = (I1 - I3) / I2;
y3 = (I2 - I1) / I3;

%torus
[u, v] = meshgrid(0:pi / 24:(2 * pi));
x = (c + a * cos(v)).* cos(u);
y = (c + a * cos(v)).* sin(u);
z = a * sin(v);

%initial rotational velocity values
wx(1) = 1;
wy(1) = 1;
wz(1) = 1;
wdotx(1) = 1;
wdoty(1) = 1;
wdotz(1) = 1;
n = 1;
%for each time step (16 seconds total)
while n < nmax
    
    %runge-kutta series based on eulers equations 
    %step 1
    kx1 = h * -y1 * wy(n) * wz(n);
    ky1 = h * -y2 * wx(n) * wz(n);
    kz1 = h * -y3 * wx(n) * wy(n);

    %step 2
    kx2 = h * -y1 * ( wy(n) + 0.5 * ky1) * ( wz(n) + 0.5 * kz1);
    ky2 = h * -y2 * ( wx(n) + 0.5 * kx1) * ( wz(n) + 0.5 * kz1);
    kz2 = h * -y3 * ( wx(n) + 0.5 * kx1) * ( wy(n) + 0.5 * ky1);

    %step 3
    kx3 = h * -y1 * ( wy(n) + 0.5 * ky2 ) * ( wz(n) + 0.5 * kz2 );
    ky3 = h * -y2 * ( wx(n) + 0.5 * kx2 ) * ( wz(n) + 0.5 * kz2 );
    kz3 = h * -y3 * ( wx(n) + 0.5 * kx2 ) * ( wy(n) + 0.5 * ky2 );

    %step 4
    kx4 = h * -y1 * ( wy(n) + ky3 ) * (wz(n) + kz3 );
    ky4 = h * -y2 * ( wx(n) + kx3 ) * (wz(n) + kz3 );
    kz4 = h * -y3 * ( wx(n) + kx3 ) * (wy(n) + ky3 );

    %Add change in velocity to previous value to get the next value
    wx(n+1) = wx(n) + kx1 / 6 + kx2 / 3 + kx3 / 3 + kx4 / 6;
    wy(n+1) = wy(n) + ky1 / 6 + ky2 / 3 + ky3 / 3 + ky4 / 6;
    wz(n+1) = wz(n) + kz1 / 6 + kz2 / 3 + kz3 / 3 + kz4 / 6;

    %set wdot values
    wdotx(n+1) = wx(n + 1);
    wdoty(n+1) = wy(n + 1);
    wdotz(n+1) = wz(n + 1);

    %
    n = n + 1;

end;

%plot 3 dimensional line to show velocity over time
plot = plot3(wdotx, wdoty, wdotz);

%set labels and title
clear title xlabel ylabel zlabel;
xlabel('\omega_x /rads^{-1}');
ylabel('\omega_y /rads^{-1}');
zlabel('\omega_z /rads^{-1}');
title('The periodic evolution of the angular speed of a torus.');
legend(plot, 'Angular Speed');
grid on;