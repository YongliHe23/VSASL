function [rf, g] = readaslprep(path)

raster.grad = 10;  % us
raster.rf = 1;     % us

maxpgwamp = 32766;

g = load([path '/grad.txt']);
mxg = 1;   % Gauss/cm
g = mxg*g/maxpgwamp;  % Gauss/cm

rho = load([path '/rho.txt']);
mxrf = 200e-3;   % Gauss
rho = mxrf*rho/maxpgwamp;

theta = load([path '/theta.txt']);
theta = pi*theta/maxpgwamp;

rf = rho.*exp(1i*theta);

return

% interpolate to 1us (rf) and 10us (gradient)
n = size(rf,1);
dt = 4;   % input raster, us
dur = n*dt;
tIn = dt/2*[(1:n)-0.5];

tOut = raster.rf*[(1:(4*n))-0.5];
rf = interp1(tIn, rf, tOut, 'linear', 'extrap');
subplot(121); plot(tOut*1e-6, abs(rf));

tOut = raster.grad*[(1:(ceil(4/10*n)))-0.5];
g = interp1(tIn, g, tOut, 'linear', 'extrap');
subplot(122); plot(tOut*1e-6, g);

