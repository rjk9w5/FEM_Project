E = 70*10^9;

p1 = [.2,	.2,	14482000000;
.5,	.5,	27667000000;
1,	1,	39234000000;
2,	2,	50642000000;
5,	5,	60705000000;
10,	10,	65017000000;
20,	20,	67417000000;
50,	50,	68944000000;
100,	100,	69456000000;
200,	200,	69733000000];

p1(:,4) = abs(p1(:,3) - E)./E *100;

figure(1)
plot(p1(1:6,2), p1(1:6,4))
title('%Error Convergence')
xlabel('L = \beta A')
ylabel('%Error of Youngs modulus')

figure(2)
plot(p1(1:6,2), p1(1:6,3)/1e9)
title('Young modulus ')
xlabel('L = \beta A')
ylabel('Youngs modulus [GPa]')


% p2 = [0.03	,8.75322;
% 0.06	,11.8588;
% 0.2	,13.1255;
% 0.6	,34.5222;
% 1	,45.5351];

p2 = [0.03	,4.97145;
0.06	,11.289;
0.2	,11.992;
0.6	,36.37;
1	,46.276];




h = log(p2(:,1));
err = log(p2(:,2));

p = diff(err)./diff(h)

sum(p)/max(size(p))


figure(3)
loglog(p2(:,1), p2(:,2))
title('Error Plot')
xlabel('ln(h)')
ylabel('ln(error)')

