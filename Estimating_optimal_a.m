clear all;
close all;
% termination criteria
err_loc = 0.;
err_glo = 0.;
tol = 1.e-5;

% geometry
a = 10;
b = 9;
d = 3.75;

ag=0;
% mesh
h = 0.5;

n=0;



%alf= 2/(1+sin(h*pi)); % overrelaxation parameter
% alf=1.9935; % for 0.00125
% alf=1.988; % for 0.0025
% alf=1.97625; % for 0.005
% alf=1.95; % for 0.01


eps_0=8.854187817e-12; % C/V/m in SI units
% boundary conditions
Vedge = 0.;
Vcentral = 100.;


alf=1.9999;
bn=10000;
x1=zeros(1,bn);
y1=zeros(1,bn);
while (n <= bn) %Can reach until alf=1.0001 thus all values within limits (1<alf<2)

% initialisation 
ist=1;
iend=round(a/h)+1;
jst=1;
jend=round(b/h)+1;
kend=round(d/h)+1;
Vold = zeros(iend,jend);
Vnew = zeros(iend,jend);
Resnew = zeros(iend,jend);
dVy = zeros(kend);
pEx= zeros(jend);
pEy= zeros(jend);
max_it = (iend*jend)*10;
% time start
cl_st=clock;

   for iter = 1:max_it  % begin iteration
    fprintf('%i\n',iter);                       
    for i = ist:iend
        for j = jst:jend
            if(j == jst && i<=kend) % capacitor plate
                Vnew(i,j) = Vcentral;
            elseif (i == iend) % external plate
                Vnew(i,j)=Vedge;
            elseif (j == jend) % external plate
                Vnew(i,j)=Vedge;
            else % inside
                if (i == ist) % vertical symmetry line 
                   vlc=Vold(i+1,j);
                else
                   vlc=Vnew(i-1,j); 
                end
                if (j == jst) % horizontal symmetry line 
                   vcd=Vold(i,j+1);
                else
                   vcd=Vnew(i,j-1);
                end
                vcu=Vold(i,j+1);
                vrc=Vold(i+1,j);
                Resnew(i,j)=0.25*(vlc+vcd+vcu+vrc)-Vold(i,j);
                Vnew(i,j)=Vold(i,j)+alf*Resnew(i,j);
            end
            tmp_loc = alf*Resnew(i,j);
            if(err_loc < tmp_loc)
                err_loc = tmp_loc;
            end
            err_glo = err_glo + abs(tmp_loc);
        end
    end
              err_glo = err_glo/(iend*jend);    

     if ( err_loc <= tol ), break, end          % check convergence
     
        Vold = Vnew;
        err_loc = 0.;
        err_glo = 0.;
   

   end

y1(1,n+1)=iter;
x1(1,n+1)=alf;
alf = alf-0.0001;
n=n+1;
fprintf('%i\n',iter);

end
plot(x1,y1);




fprintf('%i\n',alf);
fprintf('%i\n',n);


%time finish
cl_end=clock;
run_t = etime(cl_end, cl_st);

% potential plot
 for i = ist:iend
     x(i)=(i-1)*h;
 end
  for j = jst:jend
     y(j)=(j-1)*h;
  end
 
% %image(x,y,Vnew');
% contourf(x,y,Vnew',20);
xlabel('Alpha','FontSize',12);
ylabel('Number of Iterations','FontSize',12);
title('Alpha vs Number of Iterations','FontSize',12);
% %contourf(transpose(Vnew),20);
 
% charge from div method eq (6)
divE=-(Vnew(ist,jst+1)-Vnew(ist,jst))/h*h/2.;
for i = ist+1:kend
    divE=divE-(Vnew(i,jst+1)-Vnew(i,jst))/h*h;
end
divE=divE-(Vnew(kend+1,jst)-Vnew(kend,jst))/h*h/2.;
charge1=divE*eps_0;

% charge from eq (5)
% derivative from eq (14)
for i = ist:kend
    dVy(i)=(-Vnew(i,jst+2)+4*Vnew(i,jst+1)-3*Vnew(i,jst))/(2*h);
end
% trapezoidal rule
charge2=0.;
for i = ist:kend-1
    charge2 = charge2 - h*(dVy(i)+dVy(i+1))/2.*eps_0;
end
% Simpson's rule
  charge5=0.;
for i = ist:2:kend-2
    charge5 = charge5 - h*(dVy(i)+4*dVy(i+1)+dVy(i+2))/3.*eps_0;
end

% charge from energy method eq (7)
energy3 = 0.;
for i = ist:iend-1
   for j = jst:jend-1
     Eyl=-(Vnew(i,j+1)-Vnew(i,j))/h;
     Eyr=-(Vnew(i+1,j+1)-Vnew(i+1,j))/h;
     Exd=-(Vnew(i+1,j)-Vnew(i,j))/h;
     Exu=-(Vnew(i+1,j+1)-Vnew(i,j+1))/h;
     energy3 = energy3 + (0.5*h^2)*eps_0*0.5*(Eyl^2+Eyr^2+Exd^2+Exu^2);

   end
end
charge3 = energy3*2./(Vcentral-Vedge);

% charge from energy method eq (7), another way
energy4 = 0.;
for i = ist:iend-1
   for j = jst:jend-1
     Ey=0.5*(-(Vnew(i,j+1)-Vnew(i,j))/h-(Vnew(i+1,j+1)-Vnew(i+1,j))/h);
     Ex=0.5*(-(Vnew(i+1,j)-Vnew(i,j))/h-(Vnew(i+1,j+1)-Vnew(i,j+1))/h);
     energy4 = energy4 + (h^2)*eps_0*0.5*(Ey^2+Ex^2);
   end
end
charge4 = energy4*2./(Vcentral-Vedge);

% field plots
   for i = ist:iend-1
     pEy(i)=-(Vnew(i,2)-Vnew(i,1))/h;
     pEx(i)=-(Vnew(i+1,1)-Vnew(i,1))/h;
   end

% figure1 = figure;
% axes1 = axes('Parent',figure1,'YGrid','on','XGrid','on','FontSize',10);
% box(axes1,'on');
% hold(axes1,'all');
% plot(x(1:iend-1),pEy(1:iend-1),'Marker','none','LineWidth',2,'Color',[0 0 1],'DisplayName','field component Ey')
% xlabel('x','FontSize',12);
% ylabel('Ey','FontSize',12);
% title('field component Ey','FontSize',12);
% %legend('show');
% 
% figure1 = figure;
% axes1 = axes('Parent',figure1,'YGrid','on','XGrid','on','FontSize',10);
% box(axes1,'on');
% hold(axes1,'all');
% plot(x(1:iend-1)+h/2,pEx(1:iend-1),'Marker','none','LineWidth',2,'Color',[0 1 0],'DisplayName','field component Ex')
% xlabel('x','FontSize',12);
% ylabel('Ex','FontSize',12);
% title('field component Ex','FontSize',12);
% %legend('show');



%print
% fprintf('%i\n',iter);
% fprintf('%i\n',iter1);
% fprintf('%i\n',run_t);
