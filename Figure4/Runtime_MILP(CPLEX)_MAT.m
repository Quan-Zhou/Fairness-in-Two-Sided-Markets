%% Para

mu1=10;
mu2=10;

lam1=0.5;
lam2=0.5;
lam3=0;

%pay=5+rand(icount,1);
pay=[1.38,0.77,2.52,1.15,1.2,3.53,1.0,2.2,0.95,1.42,1.7,2.77,0.57,0.81,21.47,1.13];

size_range=5:14;
repeat=1:10;

%% F1 formulation inter3
run=zeros(size(size_range,2),2);
tick=1;

for t = size_range
    
icount=t;
jcount=t;
scount=2;

male=(1:2:jcount);
female=(2:2:jcount);


for r =repeat
    
dis=rand(icount,jcount);



M=binvar(icount,jcount,'full');
u=sdpvar(jcount,1);

intra5=sum((1-(u+0.01)/1.01).^2);

muf=sum((u(female)+0.01)/1.01)/size(female,2);
mum=sum((u(male)+0.01)/1.01)/size(male,2);
inter3=(mum-muf).^2;
%inter0=0;

Cust=sum(sum(M.*dis));

obj=lam1*intra5+lam2*inter0+lam3*Cust;

con=[];

for j = 1:jcount
con=[con,sum((pay(1:icount)'-dis(1:icount,j)).*M(1:icount,j))==u(j)];
con = [con,1-sum(M(:,j))>=0]; % line 4
end

for i = 1 : icount
    con=[con,sum(M(i,:))==1]; % line 3
end

ops = sdpsettings;
ops = sdpsettings('verbose',2,'solver','cplex','debug', 1);
%[model,recoverymodel] = export(con,obj,ops);

output = optimize(con,obj,ops);

run(tick,:)=[t,output.solvertime];
tick=tick+1;

end

end

%writematrix(run,'run_CPLEX_F1inter3.csv')

%% F1 formulation inter0

run=zeros(size(size_range,2),2);
tick=1;

for t = size_range
    
icount=t;
jcount=t;
scount=2;

male=(1:2:jcount);
female=(2:2:jcount);

for r =repeat
    
dis=rand(icount,jcount);

M=sdpvar(icount,jcount,'full');
u=sdpvar(jcount,1);

intra5=sum((1-(u+0.01)/1.01).^2);

%muf=sum((u(female)+0.01)/1.01)/size(female,2);
%mum=sum((u(male)+0.01)/1.01)/size(male,2);
%inter3=(mum-muf).^2;
inter0=0;

Cust=sum(sum(M.*dis));

obj=lam1*intra5+lam2*inter0+lam3*Cust;

con=[];

for j = 1:jcount
con=[con,sum((pay(1:icount)'-dis(1:icount,j)).*M(1:icount,j))==u(j)];
con = [con,1-sum(M(:,j))>=0]; % line 4
end

for i = 1 : icount
    con=[con,sum(M(i,:))==1]; % line 3
end

ops = sdpsettings;
ops = sdpsettings('verbose',2,'solver','cplex','debug', 1);
%[model,recoverymodel] = export(con,obj,ops);

output = optimize(con,obj,ops);

run(tick,:)=[t,output.solvertime];
tick=tick+1;

end

end

%writematrix(run,'run_CPLEX_F1inter0.csv')

%% Formulation 2 inter3

run=zeros(size(size_range,2),2);
tick=1;

for t = size_range
    
icount=t;
jcount=t;
scount=2;

male=(1:2:jcount);
female=(2:2:jcount);

for r =repeat
    
dis=rand(icount,jcount);

M=sdpvar(icount,jcount,'full');
u=sdpvar(jcount,1);

intra5=sum((1-(u+0.01)/1.01).^2);

muf=sum((u(female)+0.01)/1.01)/size(female,2);
mum=sum((u(male)+0.01)/1.01)/size(male,2);
inter3=(mum-muf).^2;
%inter0=0;

Cust=sum(sum(M.*dis));

obj4=sum(mu2*(sum(M,2)-1).^2,1); % mu2
%obj5=sum(mu1*(sum(M.*(pay(1:icount)'-dis),1)-u').^2,2); % mu1

obj5=0;

con=[];

for j = 1:jcount
%con=[con,sum((pay(1:icount)'-dis(1:icount,j)).*M(1:icount,j))==u(j)];
con = [con,1-sum(M(:,j))>=0]; % line 4
obj5=obj5+mu1*(sum(M(:,j).*(pay(1:icount)'-dis(:,j)))-u(j)).^2;
end

for i = 1 : icount
con=[con,sum(M(i,:))==1]; % line 3
end

obj=lam1*intra5+lam2*inter3+lam3*Cust+obj4+obj5;

ops = sdpsettings;
ops = sdpsettings('verbose',2,'solver','cplex','debug', 1);
%[model,recoverymodel] = export(con,obj,ops);

output = optimize(con,obj,ops);
run(tick,:)=[t,output.solvertime];
tick=tick+1;

end

end

writematrix(run,'run_CPLEX_F2inter3.csv')

%% Formulation 2 inter0

run=zeros(size(size_range,2),2);
tick=1;

for t = size_range
    
icount=t;
jcount=t;
scount=2;

male=(1:2:jcount);
female=(2:2:jcount);

for r =repeat
    
dis=rand(icount,jcount);

M=sdpvar(icount,jcount,'full');
u=sdpvar(jcount,1);

intra5=sum((1-(u+0.01)/1.01).^2);

%muf=sum((u(female)+0.01)/1.01)/size(female,2);
%mum=sum((u(male)+0.01)/1.01)/size(male,2);
%inter3=(mum-muf).^2;
inter0=0;

Cust=sum(sum(M.*dis));

obj4=sum(mu2*(sum(M,2)-1).^2,1); % mu2
%obj5=sum(mu1*(sum(M.*(pay(1:icount)'-dis),1)-u').^2,2); % mu1

obj5=0;

con=[];

for j = 1:jcount
%con=[con,sum((pay(1:icount)'-dis(1:icount,j)).*M(1:icount,j))==u(j)];
con = [con,1-sum(M(:,j))>=0]; % line 4
obj5=obj5+mu1*(sum(M(:,j).*(pay(1:icount)'-dis(:,j)))-u(j)).^2;
end

for i = 1 : icount
con=[con,sum(M(i,:))==1]; % line 3
end

obj=lam1*intra5+lam2*inter0+lam3*Cust+obj4+obj5;

ops = sdpsettings;
ops = sdpsettings('verbose',2,'solver','cplex','debug', 1);
%[model,recoverymodel] = export(con,obj,ops);

output = optimize(con,obj,ops);
run(tick,:)=[t,output.solvertime];
tick=tick+1;

end

end

writematrix(run,'run_CPLEX_F2inter0.csv')