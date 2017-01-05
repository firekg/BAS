load dist.txt

x=dist(:,1);
y=dist(:,2)
z=x.*y;
plot(x,z,'o')

axis([0 20 0 1])
b=polyfit(x,z,1)

2.05-b(2)-10*b(1)
% for amp: 0.03*R + 0.23 + 1.51

load dir.txt

x=dir(:,1);
y=dir(:,2)
z=x.*y;
plot(x,z,'o')

axis([0 20 0 1])
b=polyfit(x,z,1)

0.85-b(2)-10*b(1)
% for dir: 0.014*R + 0.085 + 0.626
