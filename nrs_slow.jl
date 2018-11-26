#=
    Solving Nurse Rostering Problem by Imperialist Competitive Algorithm
    Created by Zong-Sheng Wang @ 2018.05.28
=#
using CPUTime
using WZS.ICA

# Nurse Rostering Problem
const NURSES = 7
const DAYS = 7
const DUTY_TYPES = 4
const DUTY_PER_DAY = [2 2 1 2] # M E N O
const DUTY_PER_NURSE = [2 2 1 2]
const DUTY_IN_CHAR = ['m', 'e', 'n', 'o']

# Initialze Parameter of Imperialist Competitive Algorithm
const COUNTRY_NUM = 200
const IMPERIALIST_NUM = 25
const COLONY_NUM = COUNTRY_NUM - IMPERIALIST_NUM
const YEARS = 150
#const Beta = 2.0
#const Gamma = pi * 0.25
const Xi = 0.1  #ξ
const ASSIMILATION_TIMES = 50


# generates correct duty per day matrix for evaluation
function genevaldaymatrix(days, duty_types)
    M = zeros(Int16, duty_types, days)
    [ M[:, i] = DUTY_PER_DAY for i=1:days ]
    return M
end

# generates correct duty per nurse matrix for evaluation
function genevalnursematrix(nurses, duty_types)
    M = zeros(Int16, nurses, duty_types)
    [ M[i, :] = DUTY_PER_NURSE for i=1:nurses ]
    return M
end

function findfirstnminindex(A::Array, n::Int)
    base = sort(A)[n]
    I = zeros(Int, n)
    i = 1;
    for (index, value) in enumerate(A)
        if i>n break end
        if value <= base
            I[i] = index
            i+=1
        end
    end
    I
end

function calcdutynum(record)
    R = zeros(Int16, DUTY_TYPES)
    #[ R[j]=R[j]+1 for c in record for j = 1:DUTY_TYPES if c == DUTY_IN_CHAR[j] ]
    for c in record
        for j = 1:DUTY_TYPES
            if c == DUTY_IN_CHAR[j]
                R[j]+=1
            end
        end
    end
    R
end

function calcnursematrix(S)
    M = zeros(Int16, NURSES, DUTY_TYPES)
    #[ M[i,:] = calcdutynum(S[i, :]) for i=1:NURSES ]
    for i=1:NURSES
        M[i,:] = calcdutynum(S[i, :])
    end
    return M
end

function calcdaymatrix(S)
    M = zeros(Int16, DUTY_TYPES, DAYS)
    for j=1:DAYS
        M[:,j] = calcdutynum(S[:, j])
    end
    return M
end

function calcfitness(solution)
    N = calcnursematrix(solution)
    D = calcdaymatrix(solution)
    #println(N),println(D)
    return sum(abs.(EVAL_NURSE-N)) + sum(abs.(EVAL_DAY-D))
end

function gensoluation()
    rng = MersenneTwister(Base.Dates.millisecond(now())+rand(1:999))
    S = Array{Char}(NURSES, DAYS)
    rest = map((x)->x*NURSES, DUTY_PER_NURSE)
    index = 1
    while sum(rest) > 0
        while true
            r = rand(rng, 1:4)
            rest[r] == 0 && continue
            S[index] = DUTY_IN_CHAR[r]
            rest[r] -= 1
            index += 1
            break
        end
    end
    S
end

function assimilation(colony)
    R = copy(colony)
    ex_end = round.(NURSES*DAYS/2)
    exchange_num = rand(1:ex_end)
    #exchange_num = 1
    #exchange_pair_indics = zeros(Int64, exchange_num*2)
    for i=1:exchange_num
        src_index = rand(1:NURSES*DAYS)
        dest_index = rand(1:NURSES*DAYS)

        R[src_index] = colony[dest_index]
        R[dest_index] = colony[src_index]
    end
    R
end


function displaysolution(S)
    N = calcnursematrix(S)
    D = calcdaymatrix(S)
    #println(N),println(D)
    println("\t MON TUE WEB THU FRI SAT SUN m e n o")
    for i=1:NURSES
        print("$i\t")
        for j=1:DAYS
            print("  $(S[i, j]) ")
        end
        for j=1:DUTY_TYPES
            print(" $(N[i, j])")
        end
        print("\n")
    end

    for i=1:DUTY_TYPES
        print("$(DUTY_IN_CHAR[i])\t")
        for j=1:DAYS
            print("  $(D[i, j]) ")
        end
        print("\n")
    end
end

#calcuates the num of colonies for each empire
function coloniesperempire(imp_fitnesses)
    max_fitness = maximum(imp_fitnesses)
    normalized_best_fitnesses = max_fitness - imp_fitnesses
    power_of_imps = normalized_best_fitnesses / sum(normalized_best_fitnesses)
    colony_num_of_each_imp::Array{Int64} = floor(power_of_imps * COLONY_NUM)
    colony_num_of_each_imp[end] += COLONY_NUM - sum(colony_num_of_each_imp)
    return colony_num_of_each_imp
end

#assigns colonies to empire
function assigncoloniestoempire( rest, colony_indics)
    colonies = zeros(Int64, rest)
    while rest>0
        l = length(colony_indics)
        if l==0 break end
        r = rand(1:l)
        colonies[rest] = colony_indics[r]
        deleteat!(colony_indics, r)
        rest-=1
    end
    return colonies
end

function createempires(imperialists, imp_fitnesses)
    empires = Array{WZS.ICA.Empire}(IMPERIALIST_NUM)

    colony_num_of_each_imp = coloniesperempire(imp_fitnesses)
    #println(colony_num_of_each_imp)
    # assign colonies to each empire
    colony_indics = deleteat!(collect(1:COUNTRY_NUM), imperialists)
    for i=1:IMPERIALIST_NUM
        # assign colonies to imperialist
        rest = colony_num_of_each_imp[i]
        colonies = assigncoloniestoempire(rest, colony_indics)

        emp = WZS.ICA.Empire(imperialists[i], colonies)
        empires[i] = emp
    end
    return empires
end


# move colony towards imperialist
function revolution(colonies, countries, fitnesses)
    n_colony = length(colonies)
    for c in 1:n_colony
        # revolution
        c_index = colonies[c]
        for j = 1:ASSIMILATION_TIMES # 随机的交换位置，直到fitness 比之前好
            new_solution = assimilation(countries[c_index])
            new_fitness = calcfitness(new_solution)
            if new_fitness < fitnesses[c_index]
                fitnesses[c_index] = new_fitness
                countries[c_index] = new_solution
                break
            end
        end
    end
end

function calcempirepowers(empire_fitnesses)
    empire_powers = Array{Float32}(IMPERIALIST_NUM)
    sum_empires_fitnesses = sum(empire_fitnesses)
    max_empire_fitness =  maximum(empire_fitnesses)
    normalized_empire_fitnesses = max_empire_fitness - empire_fitnesses
    empire_powers = normalized_empire_fitnesses/sum_empires_fitnesses
end

function eliminating(empires, imperialists, empire_fitnesses, empire_powers)
    eliminating_empires = Array{Int64}(0)
    for i=1:length(empires)
        n_colony = length(empires[i].colonies_index)
        if n_colony == 0
            push!(empires[indmax(empire_powers)].colonies_index, empires[i].imp_index)
            push!(eliminating_empires, i)
        end
    end
    deleteat!(empires, eliminating_empires)
    deleteat!(imperialists, eliminating_empires)
    deleteat!(empire_fitnesses, eliminating_empires)
    deleteat!(empire_powers, eliminating_empires)
end

function main()
    println("Initial Parameter for ICA:
        Country Number: $COUNTRY_NUM
        Empire Number: $IMPERIALIST_NUM
        Colony Number: $COLONY_NUM
        Years: $YEARS
        Assimilations: $ASSIMILATION_TIMES\n")

    #generates inital solutions
    countries = Array{Any}(COUNTRY_NUM)
    fitnesses = zeros(Int64, COUNTRY_NUM)
    for i=1:COUNTRY_NUM
        countries[i] = gensoluation()
        fitnesses[i] = calcfitness(countries[i])
    end

    # inital Imperialists
    imperialists = findfirstnminindex(fitnesses, IMPERIALIST_NUM)
    empires = createempires(imperialists, fitnesses[imperialists[:]])
    empire_fitnesses = Array{Float32}(IMPERIALIST_NUM)
    sum_empires_fitnesses::Float32 = 0

    # begin
    for y = 1:YEARS
        for i=1:length(empires)
            colonies = empires[i].colonies_index
            n_colony = length(colonies)
            if n_colony > 0
                # move colony towards imperialist
                revolution(colonies, countries, fitnesses)
                # if new colony fitness is better than imperialist, then exchange
                #intraempirewar
                imp_index = empires[i].imp_index
                best_colony_fitness = minimum(fitnesses[colonies[:]])
                if  best_colony_fitness < fitnesses[imp_index]
                    best_index = indmin(fitnesses[colonies[:]])
                    imperialists[i] = empires[i].imp_index = colonies[best_index]
                    colonies[best_index] = imp_index
                end
            end

            #calcuates total power of empire
            empire_fitnesses[i] = fitnesses[empires[i].imp_index]
                + ((n_colony == 0) ? 0 : Xi*sum(fitnesses[colonies[:]])/n_colony)
        end

        # calcuates total normalized fitness and total power of empire
        empire_powers = calcempirepowers(empire_fitnesses)

        # eliminating empire with no colony
        eliminating(empires, imperialists, empire_fitnesses, empire_powers)

        # interempirewar
        weakest_empire_id = indmin(empire_powers)
        strongest_empire_id = indmax(empire_powers)
        weakest_colony_index = indmax(fitnesses[empires[weakest_empire_id].colonies_index[:]])
        push!(empires[strongest_empire_id].colonies_index, empires[weakest_empire_id].colonies_index[weakest_colony_index])
        deleteat!(empires[weakest_empire_id].colonies_index, weakest_colony_index)

        best_fitness = fitnesses[empires[strongest_empire_id].imp_index]
        println("YEAR: $y, best fitness: $best_fitness, rest empires: $(length(empires))")

        # STOP
        if length(empires) == 1 || y === YEARS || best_fitness == 0
            displaysolution(countries[empires[strongest_empire_id].imp_index])
            break;
        end

    end # end for yeas
end


EVAL_NURSE = genevalnursematrix(NURSES, DUTY_TYPES)
EVAL_DAY = genevaldaymatrix(DAYS, DUTY_TYPES)
@time main()
