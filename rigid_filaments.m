

clear all;


%RF = rod_functions;



%name of output file
fileLabel = 'test';

%choose shape of ribbon (1 or 2)
load('rodType_2.mat','X');
ds = mean(vecnorm(X(:,2:end) - X(:,1:end-1)));

%other parameters
h = 0.3; %radius
worms = 10; %number of filaments
diskR = 8; %radius of container
tubeL = 5; %height of the container

W0 = randCurve(X, h, worms, diskR); %initial configuration of ribbons
nw = size(W0,2)-2; %size of ribbons


%initialize dynamics

fStrength = 2; %inverse friction parameter
adhStrength = 20; %strength of adhesive force
extTorque = 20*[0,0,1]'; %external torque

%active brownian driving force parameters
driftForce = 2; 
driftNoise = 2;
zDrift = 0.1;

%contact moduli
bdryForceStrength = 5*10^2; %sets force from boundary
Kbulk = 2*10^2; %sets ribbon-ribbon contact forces






%Other parameters which can be changed
N=20000; %number of timesteps, 100000
framenumber=200; %number of frames to capture, 1000
dt= 10^(-4); %size of timesteps


%initialize contacts
V0comp = zeros(3, worms*(nw+2));
for ii = 1:worms
    V0comp(:, (ii-1)*(nw+2)+1: ii*(nw+2)) = W0(3*(ii-1)+1:3*ii,:);
end
Rw = dists(V0comp);
for ii = 1:worms-1 %get rid of ghost links between worm i head and worm i+1 tail
    Rw(:,ii*(nw+2)) = 100*h; Rw(ii*(nw+2),:) = 100*h;
end
%effective radius
iMw0 = interactMat(Rw,1.25*h);


%Initialize filament position array, orientation, force and drift direction
W = zeros(3*worms,nw+2,framenumber+1); W(:,:,1) = W0;

V0 = W0; 
qRot = quaternion(zeros(worms, 4)) + 1;
driftTh = zeros(1,worms);
Ftot = zeros(3, worms);

% main loop
tic;
for i=1:N 


    %initialize variables, timestepping takes V00 -> V0, etc.
    V00 = V0;
    for ii = 1:worms
        w00 = V00(3*ii-2:3*ii,:);
        w = rotatepoint( qRot(ii,:), w00' - mean(w00'));
        V0(3*ii-2:3*ii,:) = w' + mean(w00,2) + dt*Ftot(:, ii) + dt*driftForce*[cos(driftTh(ii)); sin(driftTh(ii)); zDrift*normrnd(0,1)];
    end


    %contact handling terms
    V0comp = zeros(3, worms*(nw+2));
    for ii = 1:worms
        V0comp(:, (ii-1)*(nw+2)+1: ii*(nw+2)) = V0(3*(ii-1)+1:3*ii,:);
    end

    Rw = dists(V0comp);
    Rw(1:nw+1,1:nw+1) = 100*h; %no self-effects for rigid filaments
    for ii = 1:worms-1 %get rid of links between ribbon i head and ribbon i+1 tail
        Rw(:,ii*(nw+2)) = 100*h; Rw(ii*(nw+2),:) = 100*h;
        
        %no self-effects for rigid filaments
        Rw( ii*(nw+2):(ii+1)*(nw+2)-1, ii*(nw+2):(ii+1)*(nw+2)-1 ) = 100*h;

    end


    %Nonlinear adhesion
    iMw1 = interactMat(Rw,1.25*h); %change in effective radius here
    iMw = (iMw0+iMw1)>0;
    pMatw = pressureThreshold(Rw,iMw,1.25*h,Kbulk,adhStrength); %change in effective radius here
    iMw0 = iMw1; %update interaction matrix for next run

    %calculate forces and torques
    Ftot = zeros(3, worms);
    iFW = iForce(V0comp,pMatw,Rw);
    for ii=1:worms

        w  =  V0(3*(ii-1)+1:3*ii,:);
        iF = iFW(:, (ii-1)*(nw+2)+1: ii*(nw+2));
        Fcont = iF;

        %random drift
        driftTh(ii) = driftTh(ii) + dt*sqrt(2*driftNoise/dt).*normrnd(0,1);

        %contact force
        Vradial = vecnorm(w(1:2,:));
        VLength = abs(w(3,:));
        dirns = w(:, Vradial>diskR-h);
        dirns = -dirns./vecnorm(dirns);

        bdryForce = 2*bdryForceStrength*dirns;
        Fcont(:, Vradial>diskR-h) = Fcont(:, Vradial>diskR-h) + bdryForce;

        dirns = w(3, VLength > tubeL-h);
        dirns = -[0*dirns; 0*dirns; sign(dirns)];

        bdryForce = 2*bdryForceStrength*dirns;
        Fcont(:, VLength>tubeL-h) = Fcont(:, VLength>tubeL-h) + bdryForce;

        Fcont = fStrength*ds*Fcont; %contact force with inverse friction parameter and length element
        Ftot(:,ii) = sum(Fcont,2); %total force on ribbon

        %contact torques
        CoR = mean(w,2);
        w_c = w - CoR; 
        iTau = sum(cross(w_c, Fcont), 2); %contact torque
        tau_tot = extTorque + iTau; %total torque

        %rotation quaternion
        th = norm(tau_tot); rotDir = tau_tot./th;
        qRot(ii,:) = quaternion([cos(dt*th/2), sin(dt*th/2)*rotDir']);

    end
    

    %save data every few timesteps
    if mod(i,N/framenumber)==0
        W(:,:,1+framenumber*i/N) = V0;

        %display timestep
        disp(i);
        
    end
end
% end main loop
LoopTime = toc;
LoopTime = LoopTime / N;

outputStruct = struct("W",W);

%save data

save(strcat(fileLabel,'.mat'), '-fromstruct', outputStruct);



%make video of simulation

%%%%%create AVI object
nFrames = size(W,3);
vidObj = VideoWriter(char(strcat(fileLabel)), 'MPEG-4');
vidObj.Quality = 100;
vidObj.FrameRate = 20;
open(vidObj);

colormap qualitative12

%%%%create movie
for j=1:2:nFrames
    clf;
    for ii = 1:worms

        %plot worm
        [x,y,z,C]=tubeplot2(W(3*(ii-1)+1:3*ii,:,j),h,20,h/5,[zeros(1,nw+2)+mod(ii,12)]);
        Aplot=surf(x,y,z,C);
        shading interp;
        set(Aplot,'meshstyle','row');
        set(Aplot, 'FaceLighting','phong','SpecularColorReflectance', 0, 'SpecularExponent', 50, 'DiffuseStrength', 1);
        
        material shiny; 

        hold on;
        
    end
    viscircles([0,0], diskR,'Color', 'k');
    axis(diskR*[-1,1,-1,1,-1,1]);

    view(2); caxis([0,11]);
    hl=camlight('left');
    
    daspect([1,1,1]); 
    grid off; box on;
    set(gca, 'XTick',[]); set(gca, 'YTick',[]); set(gca,'ZTick',[]);
    set(gca,'Color',[0.6 0.6 0.6]);
    set(gcf, 'Color', [1,1,1]);
    set(gca,'LooseInset',get(gca,'TightInset'));

    writeVideo(vidObj, getframe(gcf));
    drawnow;
    hold off;
    j

end
close(vidObj);
hold off;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function W0 = randCurve(X, h, worms, diskR)

    W0 = zeros(3*worms, size(X,2));
    for ii = 1:worms
        c = 1;
        while c>0
            Y = X;
            Y = rotatepoint(randrot, Y');
            Y = Y' + 2*diskR*[rand(2,1)-0.5; 0];
            W0(3*ii-2:3*ii,:) = Y;
            c = inContact(W0(1:3*ii,:), h);
            c = c + sum(vecnorm(Y(1:2,:)) >diskR);
        end
    end

end


function c = inContact(W0, h)

    worms = size(W0,1)/3;
    nw = size(W0,2)-2;
    %contact handling terms
    V0comp = zeros(3, worms*(nw+2));
    for ii = 1:worms
        V0comp(:, (ii-1)*(nw+2)+1: ii*(nw+2)) = W0(3*(ii-1)+1:3*ii,:);
    end

    Rw = dists(V0comp);

    Rw(1:nw+1,1:nw+1) = 100*h;
    for ii = 1:worms-1 %get rid of links between worm i head and worm i+1 tail
        Rw(:,ii*(nw+2)) = 100*h; Rw(ii*(nw+2),:) = 100*h;
        
        Rw( ii*(nw+2):(ii+1)*(nw+2)-1, ii*(nw+2):(ii+1)*(nw+2)-1 ) = 100*h;

    end

    %change in effective radius here
    iMw = interactMat(Rw,2*1.25*h);

    c = nnz(iMw);

end





function R = dists(X) %distance between i'th and j'th links
    Xm = 0.5*(X(:,1:end-1)+X(:,2:end)); %n+1 points, midpoint of each link
    L2 = diag(Xm'*Xm);
    R = -2*(Xm'*Xm) + (L2 + L2');
    R = R.^0.5; 
end





function iM = interactMat(R, h) %Set R(i,j)=0 for |i-j|<=10
    n = size(R,1)-1;
    R([1:n+2:end, n+1+1:n+2:end-1]) = 2*h; 
    R([2:n+2:end-(n+1)]) = 2*h; 
    iM = R<2*h; %iM(i,j) = 1 iff link i and j are interacting
end




function pMat = pressureThreshold(R,iM,h,Kbulk,adhStrength)
    R(iM==0) = 0;
    pMat = Kbulk*( 1*(1 - R./(2*h)) + 100*(1-R./(2*h)).^3 + 100*(1-R./(2*h)).^5 + 100*(1-R./(2*h)).^7); %set cubic term to around 100
    pMat(iM==0)=0;
    pMat(pMat<0) = -adhStrength;
end




function Fi = iForce(X, pMat, R)
    Xm = 0.5*(X(:,1:end-1)+X(:,2:end)); %n+1 points, midpoint of each link
    R(R==0) = 1;
    pMat = pMat./R;
    pMat = sum(pMat).*Xm - (Xm*pMat);
    pMat=[pMat(:,1), 0.5*(pMat(:,2:end)+pMat(:,1:end-1)), pMat(:,end)];
    Fi=pMat;
end




function [x,y,z,C]=tubeplot2(curve,r,n,ct,S)
% Usage: same as above but this gives colours the rod according to a 
% scalar field S along the curve.

  if nargin<3 || isempty(n), n=8;
     if nargin<2, error('Give at least curve and radius');
     end
  end
  if size(curve,1)~=3
    error('Malformed curve: should be [3,N]');
  end
  if nargin<4 || isempty(ct)
    ct=0.5*r;
  end

  
  %Collapse points within 0.5 r of each other
  npoints=1;
  for k=2:(size(curve,2)-1)
    if norm(curve(:,k)-curve(:,npoints))>ct
      npoints=npoints+1;
      curve(:,npoints)=curve(:,k);
    end
  end
  %Always include endpoint
  if norm(curve(:,end)-curve(:,npoints))>0
    npoints=npoints+1;
    curve(:,npoints)=curve(:,end);
  end

  %deltavecs: average for internal points.
  %           first strecth for endpoitns.
  dv=curve(:,[2:end,end])-curve(:,[1,1:end-1]);

  %make nvec not parallel to dv(:,1)
  nvec=zeros(3,1);
  [~,idx]=min(abs(dv(:,1))); nvec(idx)=1;

  xyz=zeros(3,n+1,npoints+2);
  Col=zeros(3,n+1,npoints+2); 

  %precalculate cos and sing factors:
  cfact=repmat(cos(linspace(0,2*pi,n+1)),[3,1]);
  sfact=repmat(sin(linspace(0,2*pi,n+1)),[3,1]);
  
  %Main loop: propagate the normal (nvec) along the tube
  for k=1:npoints
    convec=cross(nvec,dv(:,k));
    convec=convec./norm(convec);
    nvec=cross(dv(:,k),convec);
    nvec=nvec./norm(nvec);
    %update xyz:
    xyz(:,:,k+1)=repmat(curve(:,k),[1,n+1])+...
        cfact.*repmat(r*nvec,[1,n+1])...
        +sfact.*repmat(r*convec,[1,n+1]);
    Col(:,:,k+1)=S(k);
  end
  %finally, cap the ends:
  xyz(:,:,1)=repmat(curve(:,1),[1,n+1]);
  xyz(:,:,end)=repmat(curve(:,end),[1,n+1]);
  %,extract results:
  x=squeeze(xyz(1,:,:));
  y=squeeze(xyz(2,:,:));
  z=squeeze(xyz(3,:,:));
  Ct=squeeze(Col(1,:,:));
  C=Ct; C(:,1) = C(:,2); C(:,end) = C(:,end-1);
  %... and plot:
  if nargout<3, surf(x,y,z,C); end
end





