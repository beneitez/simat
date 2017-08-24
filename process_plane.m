% close all
clear all
fclose('all');
clc

% Reading the plane data after planes have been generated
% can also read unformatted data
% Can generate a gif of the plane
% 20160620 - Miguel Beneitez - beneitez@kth.se


plname = 'u_z0.pl';       % The name of the plane
gif    = 0;               % If we want a gif
shift  = 1;               % If the box is shifted so x0=0
laminar= 1;               % If we want to see a laminar profile in xy plane
fig    = 1;               % If we want to get some figures
interv = 20;              % Every how many intervals we want the figures
tmax   = 10000000;
filename = 'boxShift.gif'; % For the gif file

% Reading the plane files. The variables name correspond mostly to the
% names in SIMSON

fid=fopen(plname,'r','ieee-le');
eor = fread(fid,1,'int32');
Re  = fread(fid,1,'float64');
eor = fread(fid,1,'int32');
xl  = fread(fid,1,'float64');
zl  = fread(fid,1,'float64');
t   = fread(fid,1,'float64');
eor   = fread(fid,1,'float64');
eor = fread(fid,2,'int32');

nx   = fread(fid,1,'int32');
ny   = fread(fid,1,'int32');
nz   = fread(fid,1,'int32');
eor  = fread(fid,3,'int32');

% Calculate the coordinates based on the read info
% Shifting the box to have it start at x=0

xF=xl/nx*(-nx/2:1:nx/2-1)';
if shift==1
     xF = xF-xF(1); 
end

yF=(cos(pi*(0:1/(ny-1):1)))';
yF=yF+yF(1);
zF=zl/nz*(-nz/2:1:nz/2-1)';

tpl   = fread(fid,1,'int32');
ivar  = fread(fid,1,'int32');
cpl   = fread(fid,1,'float64');
fltype  = fread(fid,1,'int32');
dstar   = fread(fid,1,'float64');
eor  = fread(fid,2,'int32');

% Rescaling the boundary layer variables

xF=xF/dstar;
yF=yF/dstar;
zF=zF/dstar;
xl=xl/dstar;
zl=zl/dstar;
Re=Re*dstar;

fprintf(['Re = ',num2str(Re,'%4.2f'),'  xl = ', num2str(xl,'%4.2f'),...
    '  zl = ', num2str(zl,'%4.2f'),'  nx = ', num2str(nx,'%3i'),...
    '  ny = ', num2str(ny,'%3i'),'  nz = ', num2str(nz,'%3i'),'\n\n'])    
   
% Start to read for the different times

t = 0;
i = 0;
if tpl == 2
      uxz  = zeros(nx,nz);
   elseif tpl == 1
      uxy  = zeros(nx,ny);
end

% The loop until it is either empty or a max time is reached

while ~isempty(t) && t/dstar<tmax
    i = i+1;
   t    = fread(fid,1,'float64');
    if isempty(t)
        break
    end
    ts(i) = t/dstar;
    eor  = fread(fid,1,'float64');
    eor  = fread(fid,2,'int32');
    
    if isempty(t) 
        break
    end
    
    fprintf(['t = ',num2str(t/dstar,'%4.4f'),'\n']);
    
    if tpl == 2
        u  = fread(fid,nx*nz,'float64');
        u  = reshape(u,[nx,nz]);
        if shift==1
            aux = u(1:nx/2,:);
            u(1:nx/2,:) = u((nx/2+1):nx,:);
            u((nx/2+1):nx,:)= aux;
        end 
        uxz(:,:,i) = u;
        eor = fread(fid,2,'int32');
    elseif tpl == 1
        u  = fread(fid,nx*ny,'float64');
        u  = reshape(u,[nx,ny]);
        if shift==1
            aux = u(1:nx/2,:);
            u(1:nx/2,:) = u((nx/2+1):nx,:);
            u((nx/2+1):nx,:)= aux;
        end 
        uxy(:,:,i) = u;
        eor = fread(fid,2,'int32');
    end
end

% The planes have been read, if we now want to make a plot

if fig
    it = i;
    n  = 0;
    for i = 1:(it-1) 
        if mod(i,interv)==0
            n=n+1;
            figure(98)
            if tpl == 2
%         contourf(xF,zF,uxz(:,:,i)',20)%,'LineStyle','none')
                surf(xF,zF,uxz(:,:,i)')
                colorbar
                xlabel('x')
                ylabel('z')
                str = sprintf('xz plane, t = %f, i = %i',ts(i),i);
                title(str)
                drawnow
            elseif tpl == 1
                contourf(xF,yF,uxy(:,:,i)',10)          
                colorbar
                caxis([-0.5 0.5])
                xlabel('x')
                ylabel('y')
                str = sprintf('xy plane, t = %f',ts(i));
                title(str)
                drawnow
                
                % If we have an xy plane, there is the posibility to 
                
                if laminar
                    u_lam = uxy(nx,:,i);
                    figure(2)
                    plot(u_lam,yF)
                    axis([-0.5 1.01 0 yF(1)]) 
                    hold on
                end
            end
        
            % If we want to make our figure a gif
            
            if gif
                h = figure(98);
                axis tight manual
                frame = getframe(h);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256); 
  
                if n == 1 
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
                else 
                imwrite(imind,cm,filename,'gif','WriteMode','append'); 
                end
            end
        end
     end   
 end 
