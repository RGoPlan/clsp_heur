using ArgParse
using DelimitedFiles
using JuMP, Gurobi

const GUROBI_ENV = Gurobi.Env()

# Se precisar, comente esse bloco de código
s = ArgParseSettings()
@add_arg_table! s begin
    "--file", "-f"
        help = "matrix file path"
        arg_type = String
        required = true
end
parsed_args = parse_args(s)

function relax_and_fix(model, k=2, kl=1, tlimit=600)
    model_rf = model
    t_elapsed = 0
    α = 1
    nper = length(model[:x][:,1])
    nprod = length(model[:x][1,:])

    for var in model[:y]
        unset_binary(var)
        #=
        if is_binary(var)
            unset_binary(var)
        elseif is_integer(var)
            unset_integer(var)
        end
        =#
    end

    for i in range(1, stop=nper, step=kl)
        set_time_limit_sec(model_rf, ( (tlimit - t_elapsed) / (abs(nper - α + 1) / kl) ))
        β = min(α + k -1, nper)

        for t in α:β
            for j in 1:nprod
                #set_integer(model[:x][t, j])
                set_binary(model[:y][t, j])
                #set_integer(model[:s][t, j])
            end
        end

        t_elapsed = @elapsed optimize!(model_rf)

        if has_values(model_rf) && termination_status(model_rf) in (MOI.FEASIBLE_POINT, MOI.OPTIMAL, MOI.TIME_LIMIT)
            for t in α:(α + kl -1)
                for j in 1:nprod
                    #fix(model[:x][t,j], value.(model[:x][t,j]))
                    fix(model[:y][t,j], value.(model[:y][t,j]))
                    #fix(model[:s][t,j], value.(model[:s][t,j]))
                end
            end
            α = α + kl
            β = min(α + k -1, nper)
        else
            error("Status do modelo deu ruim")
        end
    end

    return objective_value(model_rf)
end

#parsed_args = Dict("file" => "data/540-1/1meu")
nInstancias = parse(Int, readlines(parsed_args["file"])[1])
instancias = readlines(parsed_args["file"])[2:2+nInstancias-1]

for i in 1:nInstancias
    f = readdlm(instancias[i])

    nprod = f[1,1]
    nper = f[1,2]

    inst = (
        nprod = f[1,1],
        nper = f[1,2],
        p = f[2,1],
        c = f[3,1],
        a = f[4:4+nprod-1, 1],
        b = f[4:4+nprod-1, 2],
        h = f[4:4+nprod-1, 3],
        q = f[4:4+nprod-1, 4],
        d = f[4+nprod:4+nprod+nper-1, 1:nprod]
    )

    model = Model(() -> Gurobi.Optimizer(GUROBI_ENV))
    set_optimizer_attributes(
        model, 
        "Threads" => 1,
        "MIPGapAbs" => 1e-8,
        "ImproveStartGap" => 1e-8
    )
    set_silent(model)

    @variables(model, begin
        x[t=1:inst.nper, j=1:inst.nprod], Int
        y[t=1:inst.nper, j=1:inst.nprod], Bin
        s[t=1:inst.nper, j=1:inst.nprod], Int
    end)

    @objective(
        model, Min,
        sum(inst.p*x[t,j] + inst.q[j]*y[t,j] + inst.h[j]*s[t, j] for t=1:inst.nper, j=1:inst.nprod)
    )

    # TODO: entender quando e se usar variações em t e j
    @constraints(model, begin
        capac_constr[t=1:inst.nper], sum(inst.a[j]*x[t,j] + inst.b[j]*y[t,j] for j=1:inst.nprod) <= inst.c
        stock1[j=1:inst.nprod], x[1, j] - s[1, j] == inst.d[1, j]
        stock[t=2:inst.nper, j=1:inst.nprod], inst.d[t, j] + s[t, j] == x[t, j] + s[t-1,j]
        setup[t=1:inst.nper, j=1:inst.nprod], x[t,j] <= y[t,j]*max((inst.c - inst.b[j]) / inst.a[j], sum(inst.d[t,j] for t=t:inst.nper))
        nnegatx[t=1:inst.nper, j=1:inst.nprod], x[t,j] >= 0
        nnegats[t=1:inst.nper, j=1:inst.nprod], s[t,j] >= 0
    end)

    #println(model)
    fo = relax_and_fix(model)

    println("FO: ", fo)
end
