close all;
figure;
%set a number of initial variables
a = 1;
c = 4;
M = 3.14;
h = 0.05; %timestep size in seconds
t = 16; %total time (though may be faster due to processor speed)
nmax = t / h;

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

%create data structure to store rotational velocity for each timestep
wx = zeros(1,nmax);
wy = zeros(1,nmax);
wz = zeros(1,nmax);
wdotx = zeros(1,nmax);
wdoty = zeros(1,nmax);
wdotz = zeros(1,nmax);

%create variables to store torus data 
[u,v] = meshgrid(0:pi / 24:(2 * pi));
x = (c + a * cos(v)) .* cos(u);
y = (c + a * cos(v)) .* sin(u);
z = a * sin(v);
sizex = size(x);

shape = surf(x, y, z);

%create a backup of data for later use
basex = shape.XData;
basey = shape.YData;
basez = shape.ZData;

%set initial velocity to 1
wx(1) = 1;
wy(1) = 1;
wz(1) = 1;
wdotx(1) = 1;
wdoty(1) = 1;
wdotz(1) = 1;

axis equal;
axis vis3d;    
set(gcf(), 'Renderer', 'OpenGL');
%labels and titles
xlabel('x /m');
ylabel('y /m');
zlabel('z /m');
title('Rotating Torus');
axis([-7 7 -7 7 -7 7]);

%start video 
animation = VideoWriter('animation.avi');
open(animation);

%for each time step (16 seconds total)
for n = 1:nmax
    
    %runge-kutta series based on eulers equations 
    %step 1
    kx1 = h * -y1 * wy(n) * wz(n);
    ky1 = h * -y2 * wx(n) * wz(n);
    kz1 = h * -y3 * wx(n) * wy(n);

    %Step 2
    kx2 = h * -y1 * ( wy(n) + 0.5 * ky1) * ( wz(n) + 0.5 * kz1);
    ky2 = h * -y2 * ( wx(n) + 0.5 * kx1) * ( wz(n) + 0.5 * kz1);
    kz2 = h * -y3 * ( wx(n) + 0.5 * kx1) * ( wy(n) + 0.5 * ky1);

    %Step 3
    kx3 = h * -y1 * ( wy(n) + 0.5 * ky2 ) * ( wz(n) + 0.5 * kz2 );
    ky3 = h * -y2 * ( wx(n) + 0.5 * kx2 ) * ( wz(n) + 0.5 * kz2 );
    kz3 = h * -y3 * ( wx(n) + 0.5 * kx2 ) * ( wy(n) + 0.5 * ky2 );

    %Step 4
    kx4 = h * -y1 * ( wy(n) + ky3 ) * (wz(n) + kz3 );
    ky4 = h * -y2 * ( wx(n) + kx3 ) * (wz(n) + kz3 );
    kz4 = h * -y3 * ( wx(n) + kx3 ) * (wy(n) + ky3 );

    %Add change in velocity to previous value to get the next value
    wx(n+1) = wx(n) + kx1 / 6 + kx2 / 3 + kx3 / 3 + kx4 / 6;
    wy(n+1) = wy(n) + ky1 / 6 + ky2 / 3 + ky3 / 3 + ky4 / 6;
    wz(n+1) = wz(n) + kz1 / 6 + kz2 / 3 + kz3 / 3 + kz4 / 6;
    
    %set wdot values
    wdotx(n+1) = wx(n);
    wdoty(n+1) = wy(n);
    wdotz(n+1) = wz(n);
    
    %calculate the magnitude of omega
    wmag = sqrt(wdotx(n+1)^2 + wdoty(n+1)^2 + wdotz(n+1)^2);
    
    %calculate theta
    theta = n * h * wmag;  
    
    %calculate alpha, beta gamma for rotation matrix
    hatw = [wdotx(n+1) wdoty(n+1) wdotz(n+1)] / wmag;
   
    alpha = hatw(1);
    beta = hatw(2);
    gamma = hatw(3);
    
    %rotation matrix
    lambda = [(alpha^2 * (1-cos(theta))) + cos(theta),                  (alpha * beta * (1 - cos(theta))) - (gamma * sin(theta)), (alpha * gamma * (1 - cos(theta))) + (beta * sin(theta)),
              (alpha * beta * (1 - cos(theta))) + (gamma * sin(theta)), (beta^2 * (1 - cos(theta))) + cos(theta),                 beta * gamma * (1 - cos(theta)) - alpha * sin(theta),
              (alpha * gamma * (1 - cos(theta))) - (beta * sin(theta)),  beta * gamma * (1 - cos(theta)) + (alpha * sin(theta)),  (gamma^2 * (1 - cos(theta))) + cos(theta)];
    
    %working copies of torus data 
    X = basex;
    Y = basey; 
    Z = basez;
    
    %reshape matrices so they can be multiplied by the rotation matrix
    X = reshape(X, 1, []);
    Y = reshape(Y, 1, []);
    Z = reshape(Z, 1, []);
    
    %multiply position matrix by rotational matrices
    xyz = [X; Y; Z;];
    xyzRxyz = lambda * xyz;

    %reshape matrices back to original form so they can be drawn using surf
    X = reshape(xyzRxyz(1,:), sizex);
    Y = reshape(xyzRxyz(2,:), sizex);
    Z = reshape(xyzRxyz(3,:), sizex);  
    
    %set the surf positions
    shape.XData = X;
    shape.YData = Y;
    shape.ZData = Z;
    
    %draw new position to screen and to video
    frame = getframe(gcf);  
    drawnow; 
    writeVideo(animation, frame);
end

close(animation);
