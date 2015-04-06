%%NE_155_HW6_Problem_3

int=0;
for h= [.01,.05,.1,.5,1]
int=int+1;
a=4; D=1; sig=0.2; S=8;
n=2*a/h;
A=zeros(n-1);
A(1,1)=2+h^2*sig/D;
A(1,2)=-1;
A(n-1,n-2)=-1;
A(n-1,n-1)=2+h^2*sig/D;

for j=2:n-2
        A(j,j)=2+h^2*sig/D;
        A(j,j+1)=-1;
        A(j,j-1)=-1;
end

b=(h^2*S).*ones(n-1,1);

phi_numeric=A\b;
x=[-4:h:4];
B=(-S/sig)/(exp(sqrt(sig/D)*a)+exp(-sqrt(sig/D)*a));
phi_an= @(x) B*(exp(sqrt(sig/D)*x)+exp(-sqrt(sig/D)*x))+S/sig;

phi_n=zeros(n+1,1);
phi_n(2:n)=phi_numeric;

phi_a=phi_an(x)';

e=(phi_a-phi_n)/phi_a;

emax(int)=max(max(e));

end

h= [.01,.05,.1,.5,1]';

hold on
plot(h,emax)
scatter(h,emax)

ylabel('Max Error', 'FontSize',15)
xlabel('Mesh Size (h)', 'FontSize',15)
title('Diffusion Equation Solution Error due to Mesh Size', 'FontSize', 15)


