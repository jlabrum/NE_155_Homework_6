%%NE_155_HW6_Problem_5

a=4; D=1; sig=0.7; S=8; vsig=.6; h=.1;
n=2*a/h;

%Step 1 intial values
k_0=1;
phi_0=ones(n-1,1)/norm(ones(n-1,1));

%Step 2 Elements of A
A=zeros(n-1);
A(1,1)=2*D/h^2 +sig;
A(1,2)=-D/h^2;
A(n-1,n-2)=-D/h^2;
A(n-1,n-1)=2*D/h^2 +sig;

for j=2:n-2
        A(j,j)=2*D/h^2 +sig;
        A(j,j+1)=-D/h^2;
        A(j,j-1)=-D/h^2;
end

%Step 3 Intial fission source
F= vsig*eye(n-1);
Q_0=F*phi_0;

%Step 4 a' solve for phi using GS
Diag=diag(diag(A));
L=tril(A)-Diag;
U=triu(A)-Diag;

b=(1/k_0)*Q_0;
phi_m=(-1)*inv(Diag+L)*U*phi_0+inv(Diag+L)*b;

%Step 4b' next fission source
Q_m=F*phi_m;
    
%Step 4c' next eigen value
k_m=k_0*sum(Q_m)/sum(Q_0);

%Step 4d' convergence
e1=abs(k_m-k_0);
e2=norm(abs(phi_m-phi_0));

%Repeat
iterations=1;
while e1 > 10^-4 | e2 > 10^-4
   iterations=iterations+1;
   k_0=k_m;
   Q_0=Q_m;
   b=1/k_m*Q_0;
   phi_0=phi_m;
   phi_m=(-1)*inv(Diag+L)*U*phi_0+inv(Diag+L)*b;
   Q_m=F*phi_m;
   k_m=k_0*sum(Q_m)/sum(Q_0);
   
   e1=abs(k_m-k_0);
   e2=norm(abs(phi_m-phi_0));
end

phi_m=phi_m/norm(phi_m);

iterations=iterations
k_m=k_m

phi_n=zeros(n+1,1);
phi_n(2:n)=phi_m;

x=[-4:.1:4];
hold on
plot(x,phi_n)
scatter(x,phi_n)

ylabel('\Phi  (x)', 'FontSize',15)
xlabel('x  (cm)', 'FontSize',15)
title('Solution to Eigenvalue form of Diffusion Equation', 'FontSize', 15)
