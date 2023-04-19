function estimates = fullMESIestimates(K,T,logspacing)
narginchk(1,3)
if nargin<3 || isempty(logspacing)
    logspacing = false;
end

%% load in peak derivative polynomial
folder = fileparts(mfilename('fullpath'));
load(fullfile(folder,'MESIeqDmufits.mat'),'pf0','pf1');

%% obtain first derivative and take nonuniform moving window, width = 0.5dec

sz = size(K,1,2);
LowerLimit = ones(sz)*1e-6;
row = sz(1); col = sz(2);
logT = reshape(log10(T),1,1,numel(T));

dlogT = (logT(3:end)+logT(1:end-2))/2;

if logspacing
    Kf = (movmean(K,3,3)+K)/2;  % partial noise reduction
    
    % central derivative
    dKdt = (Kf(:,:,3:end)-Kf(:,:,1:end-2))./(logT(3:end)-logT(1:end-2));
    
    % adjust initial point
    Kf = Kf(:,:,1:end-2)+(dlogT-logT(1:end-2)).*dKdt;
    
    dKdtf = min(movmean(dKdt,3,3),0);
else
    Kf = (numovmean(logT,K,logT,1,3,'tri')+K)/2; % partial noise reduction
    
    dKdt = (Kf(:,:,3:end)-Kf(:,:,1:end-2))./(logT(3:end)-logT(1:end-2));
    Kf = Kf(:,:,1:end-2)+(dlogT-logT(1:end-2)).*dKdt;
    
    dKdtf = min(numovmean(dlogT,dKdt,dlogT,1,3,'tri'),0);
end

%% estimate peak differential using second derivative

% use interpolation of second derivative
val = 1; % spacing for interpolation
[~,d2Kt] = min(dKdtf,[],3);
d2Kt(d2Kt<=val) = d2Kt(d2Kt<=val)+val;
d2Kt(d2Kt>=numel(dlogT)-val+1) = d2Kt(d2Kt>=numel(dlogT)-val+1)-val;
[rows,cols] = ndgrid(1:row,1:col);
zidx_1 = sub2ind(size(dKdtf),rows,cols,d2Kt-val);
zidx = sub2ind(size(dKdtf),rows,cols,d2Kt);
zidx1 = sub2ind(size(dKdtf),rows,cols,d2Kt+val);

b1 = dKdtf(zidx_1);
b2 = (dKdtf(zidx)-dKdtf(zidx_1))./(dlogT(d2Kt)-dlogT(d2Kt-val));
b3 = ((dKdtf(zidx1)-dKdtf(zidx))./(dlogT(d2Kt+val)-dlogT(d2Kt))-b2)./(dlogT(d2Kt+val)-dlogT(d2Kt-val));

logTmax = min(max((b3.*(dlogT(d2Kt-val)+dlogT(d2Kt))-b2)./b3/2,dlogT(d2Kt-val)-1),dlogT(d2Kt+val)+1);

Tmax = 10.^logTmax;
dKmax = b1+b2.*(logTmax-dlogT(d2Kt-val))+b3.*(logTmax-dlogT(d2Kt-val)).*(logTmax-dlogT(d2Kt));

%% estimate K0 (beta) and Kinf (beta*(1-rho)^2)

% estimate K0 (beta)
idx = 1;
x = 10.^dlogT(idx)./Tmax;
srx = sqrt(x);
sr3x = x.^(1.5);
x2 = x.^2;
B_K0 = -sqrt(1-0.045106*srx+1.8508*x-1.1288*sr3x+1.3459*x2)/log(10);
B_K1 = -sqrt(1+1.03401*srx+0.997606*x+0.08485*sr3x+0.0942633*x2)*2/log(10);
K0 = min(max(Kf(:,:,idx)+B_K0.*dKdtf(:,:,idx),0),1);

% estimate beta*(1-rho)^2
idx = numel(dlogT);
x = 10.^dlogT(idx)./Tmax;
srx = sqrt(x);
sr3x = x.^1.5;
K_P0 = sqrt(1-0.00128996*srx+1.16603*x-0.17042*sr3x+(0.9182*x).^2)./(0.9182*x)/log(10);
K_P1 = sqrt(sqrt(1+0.98233*srx+0.754227*x+0.04413*sr3x+(0.23324*x).^2)./(0.23324*x))/log(10);
Kinf0 = min(max(Kf(:,:,idx)+K_P0.*dKdtf(:,:,idx),0),1);
rho0 = (1-sqrt(Kinf0./K0));

% find peak slope location between theoretical ranges
Dm00 = K0.*polyval(pf0,rho0); % fit Dmu = 0, estimates assume Dmu = 0
Dm10 = K0.*polyval(pf1,rho0); % fit Dmu = 1, estimates assume Dmu = 0
dm = (min(max((dKmax-Dm00)./(Dm10-Dm00),0),1));%>(dKmax-Dm11)./(Dm01-Dm11);
% figure
% histogram(dm)
dm = (sqrt(2)/2*abs(dm-0.5).^0.5.*sign(dm-0.5)+0.5);
% figure
% histogram(dm)

%% re-estimate beta and rho using known dm

% estimate K0 (beta)
B_K = 1./(dm./B_K1+(1-dm)./B_K0);
K0 = min(max(Kf(:,:,1)+B_K.*dKdtf(:,:,1),0),1);

% estimate beta*(1-rho)^2
K_P = 1./(dm./K_P1+(1-dm)./K_P0);
Kinf = min(max(Kf(:,:,end)+K_P.*dKdtf(:,:,end),0),1);
rho = (1-sqrt(max(Kinf./K0,0)));

% re-evaluate beta with max slope for underestimates again
beta = K0;

%% use determined beta and rho for Kcrit estimation (T==tc)

Kcrit = min(beta.*(rho.^2.*(0.5-dm+(6*dm+0.5)/exp(2))+...
    rho.*(1-rho).*((108*dm+4)/exp(1)-40*dm)+(1-rho).^2),0.8*beta);

[~,Kclose] = max((K>=Kcrit)./-logT,[],3);
LowerLimit(Kclose==numel(logT)) = T(end);
Kclose(Kclose==numel(logT)) = numel(logT)-1;
zidx = sub2ind(size(K),rows,cols,Kclose);
zidx1 = sub2ind(size(K),rows,cols,Kclose+1);

tc = (Kcrit-K(zidx)).*(logT(Kclose+1)-logT(Kclose))./(K(zidx1)-K(zidx))+logT(Kclose);

tc = min(max((10.^tc),LowerLimit),1);

estimates(:,:,1) = tc;
estimates(:,:,2) = beta;
estimates(:,:,3) = rho;
estimates(:,:,4) = dm;
%% end of initial code

%% correction for rho
% Ke = Pipeline.fullMESIeq(cat(3,tc,beta,rho,dm,zeros(row,col)),10.^dlogT([1 end]));
% c_offset = Kf(:,:,end)-Ke(:,:,end);
% Kinf = Kinf+c_offset;
%
% beta = K0;
% rho2 = (1-sqrt(max(Kinf./K0,0)));
% 
% Kcrit = min(beta.*(rho2.^2.*(0.5-dm+(6*dm+0.5)/exp(2))+...
%     rho2.*(1-rho2).*((108*dm+4)/exp(1)-40*dm)+(1-rho2).^2),0.8*beta);
% 
% 
% Kdist = Kf-Kcrit;
% [~,idx] = max((Kf>Kcrit).*logT,[],3);
% [~,Kclose] = min(abs(Kdist),[],3);
% 
% Kclose(Kclose==numel(logT)) = numel(logT)-1;
% zidx = sub2ind(size(K),rows,cols,Kclose);
% zidx1 = sub2ind(size(K),rows,cols,Kclose+1);
% 
% tc2 = (Kcrit-K(zidx)).*(logT(Kclose+1)-logT(Kclose))./(K(zidx1)-K(zidx))+logT(Kclose);
% tc2 = min(max((10.^tc2),1e-8),1);