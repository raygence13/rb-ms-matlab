function m3=mm(n)
%
m2= zeros(n,n+1);
m3= zeros(n,n+1,n+3);

% make column vector length n
m1= (1:n)';
m1

% add columns at right to make matrix
for m=1:(n+1),
 m2(:,m)+= m1+2*(m-1);
 end;
m2

% add matrixes at right to make 3d matrix
for m=1:(n+3),
 m3(:,:,m)+= m2+(n*(n+1))*(m-1);
 end; 
size(m3)

end;

