%%NE_155_HW6_Problem_4a
int=0;
iterations_J=0;
iterations=zeros(8,1);
for error= [10^-3, 10^-5]
for h=[1,.5,.1,.05]
    
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

x_0=zeros(n-1,1);

Diag=diag(diag(A));

x_k=inv(Diag)*(Diag-A)*x_0+inv(Diag)*b;

e=norm(abs(x_k-x_0));

while e > error
    x_k=inv(Diag)*(Diag-A)*x_0+inv(Diag)*b;
    e=norm(abs(x_k-x_0));
    x_0=x_k;
    iterations_J=iterations_J+1;   
end
iterations(int)=iterations_J;

end
end

h= [1,.5,.1,.05]';

hold on
plot(h,iterations(1:4),'r')
plot(h,iterations(5:8))

ylabel('Iterations', 'FontSize',15)
xlabel('Mesh Size (h)', 'FontSize',15)
title('Diffusion Equation Using Jacobi Method', 'FontSize', 15)
legend('Error=10^-3','Error=10^-5')


scatter(h,iterations(1:4),'g')
scatter(h,iterations(5:8),'g')