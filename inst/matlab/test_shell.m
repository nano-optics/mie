clearvars
close all

nNmax=50;

lambda=transpose(300:1:600);  

Ca={30,31};
Cepsilon={epsilon_Ag(lambda),1.34^2,1.33.^2};


stMieRes=mie_FF(lambda,Cepsilon,Ca,nNmax);

Qext = stMieRes.Qext; 
Qsca = stMieRes.Qsca; 
Qabs = stMieRes.Qabs; 



plot(lambda, [Qext, Qsca, Qabs])