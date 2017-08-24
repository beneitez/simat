% clear all
% clc
%
%   Script for reading a pressure plane from the DNS
%   Miguel Beneitez - beneitez@kth.se 16082016
%

function [p,xF,yF,zF,Lx,Ly,Lz,t,Re,fltype,dstar]=pread(filename)

fid=fopen(filename,'r','ieee-le');
eor = fread(fid,1,'int32');
Re  = fread(fid,1,'float64');
eor = fread(fid,1,'int32');
Lx  = fread(fid,1,'float64');
Lz  = fread(fid,1,'float64');
t   = fread(fid,1,'float64');
eor = fread(fid,1,'float64');
eor = fread(fid,2,'int32');
NNx  = fread(fid,1,'int32')/2;
NNy  = fread(fid,1,'int32');
NNz  = fread(fid,1,'int32');
nsyn = fread(fid,1,'int32');
eor = fread(fid,2,'int32');
fltype = fread(fid,1,'int32');
dstar  = fread(fid,1,'float64');
eor    = fread(fid,1,'int32');

Ly = 2/dstar;

realpos=kron(ones(1,NNx),[1 0]);


disp(' - Reading p');
for indz=1:NNz
  for indy=1:NNy
    fread(fid,1,'int');
    vec=fread(fid,NNx*2,'float64');
    rlu(:,indz,indy)=vec(~~realpos);
    ilu(:,indz,indy)=vec(~realpos);
    fread(fid,1,'int');
  end
end

fclose(fid);

scale=1/dstar;
padx=0;
padz=0;
NxF=2*NNx;  
NzF=NNz;
xF=Lx/NxF*(-NxF/2:1:NxF/2)';
zF=Lz/NzF*(-NzF/2:1:NNz/2)';
yF=scale*(1+cos(pi*(0:1/(NNy-1):1)))';
xF=-xF(1)+xF;

%
% Shift velocity field in the streamwise direction in order
% to move the fringe to the end of the domain for spatial
% flows
%
kxvec=linspace(0,2*pi/Lx*(NNx-1),NNx);
kzvec=linspace(0,2*pi/Lz*(NNz/2-1),NNz/2);
kzvec=[kzvec -fliplr(kzvec(1:end))-1.];

xs = Lx/2.;
zs = 0.;

for i=1:NNx
  argx = -xs*kxvec(i);
  cx(i) = cos(argx);
  sx(i) = sin(argx);
end
for k=1:NNz
  argz = -zs*kzvec(k);
  for i=1:NNx
    ca(i)=cx(i)*cos(argz)-sx(i)*sin(argz);
    sa(i)=cx(i)*sin(argz)+sx(i)*cos(argz);
  end
  for j=1:NNy
    for i=1:NNx
      hr=rlu(i,k,j)*ca(i)-ilu(i,k,j)*sa(i);
      ilu(i,k,j)=ilu(i,k,j)*ca(i)+rlu(i,k,j)*sa(i);
      rlu(i,k,j)=hr;
    end
  end
end

p=reshape(complex(rlu,ilu),NNx,NNz,NNy);
vel=ccat(3,p,p,p);
vel=ccat(2,vel(:,1:NNz/2,:),vel(:,NNz/2+2:end,:));
[phys,NNx,NNy,NNz]=fou2phys(vel,0,0);   
p=phys(:,:,1+0*NNy:1*NNy);
