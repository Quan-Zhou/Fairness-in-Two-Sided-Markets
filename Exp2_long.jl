#using Pkg
#Pkg.add("CSV")
#Pkg.add(path="https://github.com/wangjie212/TSSOS")

using TSSOS
using DynamicPolynomials
using CSV
using DataFrames
using Statistics
using Random

trip_distance=DataFrame(CSV.File("taxi2000.csv"))[!,"trip_distance"]
#num_sample=2000
#trip_distance=CSV.File("taxi2000.csv",header=1,limit=num_sample,threaded=false)[2000];

# new formulation 3 07/11/22
# calculate 3 objective: intra-fair inter-fair customer-careF3_tssos(5,5,taxi[1:5],0.5,0.5,0)
function F3_tssos_21(icount,jcount,taxi,lam1,lam2,lam3)
    #icount=5 # the number of customers
    #jcount=5 # the number of drivers
    #female=[2,4]
    #male=[1,3,5]
    female=collect(2:2:jcount)
    male=collect(1:2:jcount)
    
    mu1=10
    mu2=10
    mu3=0

    #lam1=1
    #lam2=1
    #lam3=0
    
    dis=rand(icount,jcount)*0.5
    
    # set nc operators
    @polyvar u[1:jcount] M[1:icount*jcount]; #z[1:jcount]; 
    var=vcat(u,M);

    #l2=[z[j]-1+(u[j]+0.01)/1.01 for j in 1:jcount];
    #l3=[z[j]+1-(u[j]+0.01)/1.01 for j in 1:jcount];
    #l4=[b1-(1/length(female))*sum(z[j] for j in female),b1-(1/length(male))*sum(z[j] for j in male)];
    #l5=[b2-sum(M[(i-1)*(jcount)+j]*dis[i,j] for j in 1:jcount) for i in 1:icount];
    l6=[1-sum(M[(i-1)*(jcount)+j] for i in 1:icount) for j in 1:jcount];
    #l7=[M[(i-1)*(jcount)+j]*(M[(i-1)*(jcount)+j]-1) for i in 1:icount for j in 1:jcount]; #equality #numeq=25
    l8=[u[j] for j in 1:jcount];
    l9=[M[(i-1)*(jcount)+j] for i in 1:icount for j in 1:jcount];

    obj=lam1*sum((1-(u[j]+0.01)/1.01)^2 for j in 1:jcount)+lam2*(sum((u[j]+0.01)/1.01 for j in female)-sum((u[j]+0.01)/1.01 for j in male))^2+lam3*sum(M[(i-1)*(jcount)+j]*dis[i,j] for j in 1:jcount for i in 1:icount)+sum(mu1*(sum(M[(i-1)*(jcount)+j]*(taxi[i]-dis[i,j]) for i in 1:icount)-u[j])^2  for j in 1:jcount)+sum(mu2*(sum(M[(i-1)*(jcount)+j] for j in 1:jcount)-1)^2  for i in 1:icount)+sum(M[(i-1)*(jcount)+j]*(M[(i-1)*(jcount)+j]-1) for i in 1:icount for j in 1:jcount)

    # pop
    pop=vcat(obj,l6); #,l7,l8,l9

    d=1; # the relaxation order
    opt,sol,data=cs_tssos_first(pop,var,d,numeq=jcount,TS="MD",solution=true);
    #opt,sol,data=cs_tssos_higher!(data,TS="MD")
    
    u_sol=sol[1:jcount]
    M_sol=sol[jcount+1:end]
    #intra5=sum((1-(u_sol[j]+0.01)/1.01)^2 for j in 1:jcount)
    
    mu=sum((u_sol[j]+0.01)/1.01 for j in 1:jcount)/jcount
    muf=sum((u_sol[j]+0.01)/1.01 for j in female)/length(female)
    mum=sum((u_sol[j]+0.01)/1.01 for j in male)/length(male)
    
    intra4=1/(2*mu*jcount^2)*sum(abs((u_sol[k]+0.01)/1.01-(u_sol[j]+0.01)/1.01) for k in 1:jcount for j in 1:jcount)
    care=sum(M_sol[(i-1)*(jcount)+j]*dis[i,j] for j in 1:jcount for i in 1:icount)
    
    intra2=abs(1/jcount*sum(((u_sol[j]+0.01)/(1.01*mu))*log((u_sol[j]+0.01)/(1.01*mu)) for j in 1:jcount))
    intra3=abs(1/jcount*sum(log((1.01*mu)/(u_sol[j]+0.01)) for j in 1:jcount))
    
    inter1=abs((length(female)/jcount)*muf/mu*log(muf/mu)+(length(male)/jcount)*mum/mu*log(mum/mu))
    inter2=abs((length(female)/jcount)*log(mu/muf)+(length(male)/jcount)*log(mu/mum))
    inter3=abs(sum((u_sol[j]+0.01)/1.01 for j in female)-sum((u_sol[j]+0.01)/1.01 for j in male))
    
    return [intra2,intra3,intra4,inter1,inter2,inter3,care]
    #return [intra4,inter3,care]
end

fold=10
# for pareto optimal
icount=10 # the number of customers
jcount=20 # the number of drivers
repeat=5
lam1_range=0.5:0.1:1

for f in 1:fold
    #lam2_range=1:-0.1:0 #step 0.1
    df_212R = DataFrame(Intra_Fair2=Float64[],Intra_Fair3=Float64[],Intra_Fair4=Float64[],Inter_Fair1=Float64[],Inter_Fair2=Float64[],Inter_Fair3=Float64[],Customer_Care=Float64[],lam1=Float64[],lam2=Float64[],lam3=Float64[])
    trip=trip_distance[randperm(size(trip_distance)[1])[1:jcount]]
    for lam1 in lam1_range
        #for lam2 in 0:0.1:(1-lam1)
        if lam1 == 1
            multi=size(lam1_range)[1]-1
        else
            multi=1
        end
        for r in 1:repeat*multi
            obj=F3_tssos_21(icount,jcount,trip,lam1,1-lam1,0) #1-lam1-lam2
            push!(df_212R, push!(obj,lam1,1-lam1,0))
        end
        #end
    end
    filename=string("ParetoOptimal_R_long",f,".csv")
    CSV.write(filename, df_212R,header=true)
end

# for pareto optimal
icount=10 # the number of customers
jcount=20 # the number of drivers
repeat=1
lam1=0.5
num_trials=100
#lam2_range=1:-0.1:0 #step 0.1
df_random = DataFrame(Intra_Fair2=Float64[],Intra_Fair3=Float64[],Intra_Fair4=Float64[],Inter_Fair1=Float64[],Inter_Fair2=Float64[],Inter_Fair3=Float64[],Customer_Care=Float64[],Trial=Int64[])

for trial in 1:num_trials
    trip=trip_distance[randperm(size(trip_distance)[1])[1:jcount]]
    for r in 1:repeat
        obj=F3_tssos_21(icount,jcount,trip,lam1,1-lam1,0) #1-lam1-lam2
        push!(df_random, push!(obj,trial))
    end
end
CSV.write("ParetoOptimal_random_long.csv", df_random,header=true)
