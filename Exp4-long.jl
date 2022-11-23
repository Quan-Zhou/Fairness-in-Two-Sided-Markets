using TSSOS
using DynamicPolynomials
using CSV
using DataFrames
using Statistics
using Random

trip_distance=DataFrame(CSV.File("taxi2000.csv"))[!,"trip_distance"]

# new formulation 3 07/11/22
# calculate 3 objective: intra-fair inter-fair customer-careF3_tssos(5,5,taxi[1:5],0.5,0.5,0)
function F3_tssos_21(icount,jcount,taxi,lam1,lam2,lam3)

    female=collect(2:2:jcount)
    male=collect(1:2:jcount)
    
    mu1=10
    mu2=10
    mu3=0
    
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
    runtime=@elapsed cs_tssos_first(pop,var,d,numeq=jcount,TS="MD",solution=false);
    return runtime
end

runInter0=F3_tssos_21(12-10,12,trip_distance[1:20],1,0,0);

# for runtime 
#push!(df, push!(obj,lam1,lam2,1-lam1-lam2))
repeat=5
size_range=12:1:25
df_run = DataFrame(size=Int64[],runInter0=Float64[],runInter3=Float64[])
for t in size_range
    for r in 1:repeat
        trip=trip_distance[randperm(size(trip_distance)[1])[1:t]]
        runInter0=F3_tssos_21(t-10,t,trip,1,0,0);
        runInter3=F3_tssos_21(t-10,t,trip,0.5,0.5,0);
        push!(df_run,(t,runInter0,runInter3))
    end
end

CSV.write("runtime_SDP_long.csv", df_run,header=true)