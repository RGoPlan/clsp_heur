using ArgParse
using DelimitedFiles
using JuMP, Gurobi

include("modules/Heuristics.jl")
using .Heuristics

include("modules/TTMcCinstance.jl")
using .Instances

# define o argumento --file no terminal
s = ArgParseSettings()
@add_arg_table! s begin
    "--file", "-f"
        help = "list file path"
        arg_type = String
        required = true
end
parsed_args = parse_args(s)


#parsed_args = Dict("file" => "data/540-1/1meu")
instancias = readlines(parsed_args["file"])

open("log", "a") do io
    write(io, "Instance & J & T & rfUB & rfIter & rfBT & rfTime & foUB & foImprov & foIter & foRnds & foTime & status & totalTime \\\\ \n")
end

function fs(n)
    x = round(float(n), digits=4)
    return x
end

for i in 1:length(instancias)

    # lê a instância
    println("$(i)/$(length(instancias)) - Rodando instância $(instancias[i])")
    inst = read_instance_TTMcC(instancias[i])

    # modela a instância
    model = model_TTMcC(inst)
    set_silent(model)
    set_optimizer_attributes(
        model,
        "Threads" => 1,
        "IntFeasTol" => 1e-6,
        "MIPGap" => 1e-8,
        "LogToConsole" => 1 #TODO: testar a diferença para set_silent
    )

    # aplica o relax-and-fix e o fix-and-optimize
    model_rf, n_inter_rf, n_bt, t_elap_rf = relax_and_fix(model)
    
    if has_values(model_rf) && termination_status(model_rf) in (MOI.FEASIBLE_POINT, MOI.OPTIMAL, MOI.TIME_LIMIT)
        obj_value_rf = objective_value(model_rf)

        model_fo, n_inter_fo, n_rounds, t_elap_fo = fix_and_optimize(model_rf)
        obj_value_fo = objective_value(model_fo)

        improv = abs(obj_value_fo - obj_value_rf) / obj_value_rf * 100
    else
        obj_value_rf, obj_value_fo, improv, n_inter_fo, n_rounds, t_elap_fo = 0, 0, 0, 0, 0, 0
    end
    
    time = t_elap_rf + t_elap_fo

    open("log", "a") do io
        write(io, "$(instancias[i]) & $(fs(length(model[:x][1,:]))) & $(fs(length(model[:x][:,1]))) & $(fs(obj_value_rf)) & $(fs(n_inter_rf)) & $(fs(n_bt)) & $(fs(t_elap_rf)) & $(fs(obj_value_fo)) & $(fs(improv)) & $(fs(n_inter_fo)) & $(fs(n_rounds)) & $(fs(t_elap_fo)) & $(termination_status(model_rf)) & $(fs(time)) \\\\ \n")
    end
end

open("log", "a") do io
    write(io, "\n\n\n")
end
