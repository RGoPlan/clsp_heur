using ArgParse
using DelimitedFiles
using JuMP, Gurobi

include("modules/Heuristics.jl")
using .Heuristics

include("modules/TTMcCinstance.jl")
using .Instances

#= define o argumento --file no terminal
s = ArgParseSettings()
@add_arg_table! s begin
    "--file", "-f"
        help = "list file path"
        arg_type = String
        required = true
end
parsed_args = parse_args(s)
=#

parsed_args = Dict("file" => "data/540-1/inputListCG")
instancias = readlines(parsed_args["file"])

for i in 1:length(instancias)
 
    # lê a instância
    inst = read_instance_TTMcC(instancias[i])

    # modela a instância
    model = model_TTMcC(inst)
    set_silent(model)
    set_optimizer_attributes(
        model, 
        "Threads" => 1,
        "MIPGapAbs" => 1e-8,
        "ImproveStartGap" => 1e-8
    )

    # aplica o relax-and-fix e o fix-and-optimize
    model_rf = relax_and_fix(model)
    model_fo = fix_and_optimize(model_rf)

    println("R&F: ", objective_value(model_rf), " | F&O: ", objective_value(model_fo))
end
