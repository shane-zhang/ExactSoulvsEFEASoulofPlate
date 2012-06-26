% 1/3 octave band, central freq. 125 250 500 1000 2000 4000 8000 16000
% 800hz
clear
clc
%This script is used to calculate the modal energy density-analtical
%solution
format long 
syms aa m n x y w1 w2 W f w_f 
m1 = 3;
n1 = 3;
F=1;
u=0.2;
v=0.2;
a=1;
b=1;
eta=0.1;
E=7.1E10;
mu=0.3;
h=0.001;
D=(E*h^3)/(12*(1-mu^2))
rho=2700;
w1=(4*F*sin(m*pi*u/a)*sin(n*pi*v/b)*sin(m*pi*x/a)*sin(n*pi*y/b))/(a*b*(1+j*eta)*pi^4*D*((m/a)^2+(n/b)^2)-rho*h*(2*pi*f)^2)
wk=0
jis1=1
for mm=1:1:m1
    w1=subs(w1,m,mm)
    wk=wk+w1
    jis1=jis1+1
    w1=(4*F*sin(m*pi*u/a)*sin(n*pi*v/b)*sin(m*pi*x/a)*sin(n*pi*y/b))/(a*b*(1+j*eta)*pi^4*D*((m/a)^2+(n/b)^2)-rho*h*(2*pi*f)^2)
end
wc=0
jis2=0
for nn=1:1:n1
    wkk=subs(wk,n,nn)
    wc=wc+wkk
    wkk=wk
    jis2=jis2+1
end
w=wc

D1=diff(w,'x',2);D2=diff(w,'y',2);D3=diff(diff(w,'x'),'y');
W = 0.25*D*((abs(D1))^2+(abs(D2))^2+2*mu*real(D1*conj(D2))+2*(1-mu)*(abs(D3))^2+rho*(2*pi*f)^2*(abs(w))^2/D) 
fc=800; %central frequency
fu=1400; %upper frequency
fl=1120; %lower frequency
cg=2*(((2*pi*fc)^2*(D/(rho*h)))^0.25)
deltaf=cg^2/(8*pi*fc*a*b);%frequency increment
N=floor(((fu-fl)/deltaf)+1)  %incrent count
N=3
format long 
xfoot = 1
jis=1
% ff=fl:deltaf:fu;
%  w_av=0;
 for x1=0:0.1:1
    yfoot = 1
    for y1=0:0.1:1
        
        ff=fl:deltaf:fu;
         w_av=0;
    for ii=1:N
    f=ff(ii);
    w_f=abs(subs(eval(W),{x,y},{x1,y1}));
    w_av=w_av+w_f;
    end
w_av_t(xfoot,yfoot)=w_av*deltaf/(fu-fl)
		yfoot = yfoot + 1
        jis=jis+1
    end
	xfoot = xfoot + 1
end
% for ii=1:N
%     f=ff(ii);
%     w_f=eval(W);
%     w_av=w_av+w_f;
% end
% w_av=w_av*deltaf/(fu-fl)

% format long 
% xfoot = 1
% jis=1
% for x1=0:0.1:1
%     yfoot = 1
%     for y1=0:0.1:1
%         WA(xfoot,yfoot)=abs(subs(w_av,{x,y},{x1,y1}));
% 		yfoot = yfoot + 1
%         jis=jis+1
%     end
% 	xfoot = xfoot + 1
% end
[xa,ya]=meshgrid(1:11);
xa=xa*0.1-0.1;
ya=ya*0.1-0.1;
mesh(xa,ya,double(w_av_t))