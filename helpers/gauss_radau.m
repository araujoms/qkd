function [w t] = gauss_radau(m)

J = zeros(m,m);
for n=1:m-1
	J(n,n) = 0.5;
	J(n,n+1) = n/(2*sqrt(4*n^2-1));
	J(n+1,n) = J(n,n+1);
end
J(m,m) = (3*m-1)/(4*m-2);

[v d] = eig(J);

w = (v(1,:).^2)';
t=diag(d);
	

