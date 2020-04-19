# data
groups = 3   # we need to finally divide them into 3 groups

table = readcsv("dairy9.csv") # read data from file
col_names = table[1,:]         # first row - ID, LAC, DIM, MY, FAT, BW

values = table[2:end, :]
number_cows = length(values[:,1]) # total number of cows in the data
println("Total number of cows = ", number_cows)

cow_id = values[:,1]  # cows unique ids
LAC = values[:,2]     # LAC values
DIM = values[:,3]     # DIM values
MY = values[:,4]      # MY values
FAT = values[:,5]     # FAT values
BW = values[:,6]      # BW values

DMI = zeros(number_cows)
NE = zeros(number_cows)
CP = zeros(number_cows)
for i=1:number_cows
    DMI[i] = (0.372*(0.4+0.15*FAT[i])*MY[i] + 0.0968*BW[i]^0.75)*(1-exp(-0.192*((DIM[i]/7)+3.67)))
    NE[i] = (0.079*BW[i]^0.75 + MY[i]*(0.36+0.0969*FAT[i]))/DMI[i]
    CP[i] = ((104.78+0.73*BW[i]-0.00015432*BW[i]^2) + (MY[i]*(4586+1036*FAT[i]))/100)*0.001/DMI[i]
end

using JuMP, AmplNLWriter, CoinOptServices

m = Model(solver=BonminNLSolver(["bonmin.nlp_log_level=0"; "bonmin.bb_log_level=0"]))
@variable(m, 0 <= cow_group[1:number_cows, 1:groups] <= 1, Int)  # matching of a cow id to a group number
@variable(m, ne_level[1:groups] >= 0)
@variable(m, cp_level[1:groups] >= 0)
for i=1:number_cows
    @constraint(m, sum(cow_group[i,:]) == 1)   # each cow should be assigned to one group only
end
for i=1:groups
    @constraint(m, sum(cow_group[:,i]) == number_cows/groups) # in each group the number of cows is 30/3=10
end

for k=1:groups
    @constraint(m, ne_level[k] == 1.3*sum(NE[i]*cow_group[i,k] for i=1:number_cows)/3)
    @constraint(m, cp_level[k] == 1.3*sum(CP[i]*cow_group[i,k] for i=1:number_cows)/3)
end

# minimize this objective
@NLobjective(m, Min, 0.07*sum(sum(DMI[i]*cow_group[i,g]*ne_level[g] for i=1:number_cows) for g=1:groups) + 
                    0.4*sum(sum(DMI[i]*cow_group[i,g]*cp_level[g] for i=1:number_cows) for g=1:groups))
status = solve(m)
                                                            
# Printing the output of the program
println("Groupings of cows:")
println("ID g-1 g-2 g-3")
for i=1:number_cows
    if cow_id[i]<=9
        print(" ", cow_id[i]," ", getvalue(cow_group[i,1]), " ", getvalue(cow_group[i,2])," ", getvalue(cow_group[i,3]),"\n")
    else 
        print(cow_id[i]," ", getvalue(cow_group[i,1]), " ", getvalue(cow_group[i,2])," ", getvalue(cow_group[i,3]),"\n") 
    end
end
println("Objective value = ", getobjectivevalue(m))
