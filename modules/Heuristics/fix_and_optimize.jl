using JuMP

function fix_and_optimize(model, k=2, kl=1, tlimit=600)
    model_fo = model
    t_elapsed = 0
    α = 1
    nper = length(model[:x][:,1])
    nprod = length(model[:x][1,:])

    for i in range(1, stop=nper, step=kl)

        # define a janela e o tempo limite de sua execução
        set_time_limit_sec(model_fo, ( (tlimit - t_elapsed) / (abs(nper - α + 1) / kl) ))
        β = min(α + k -1, nper)
        model_temp = model_fo

        # retira os valores fixos das variáveis y da janela
        for t in α:β
            for j in 1:nprod
                #unfix(model[:x][t,j])
                unfix(model_temp[:y][t,j])
                #unfix(model[:s][t,j])
            end
        end

        # seta as variáveis y da janela como binária
        for t in α:β
            for j in 1:nprod
                #set_integer(model[:x][t, j])
                set_binary(model_temp[:y][t, j])
                #set_integer(model[:s][t, j])
            end
        end

        # resolve o problema e armazena o tempo de execução
        t_elapsed = @elapsed optimize!(model_temp)

        for t in α:β
            for j in 1:nprod
                #fix(model[:x][t,j], value.(model[:x][t,j]))
                fix(model_temp[:y][t,j], value.(model_temp[:y][t,j]))
                #fix(model[:s][t,j], value.(model[:s][t,j]))
            end
        end

        if has_values(model_temp) && termination_status(model_temp) in (MOI.FEASIBLE_POINT, MOI.OPTIMAL, MOI.TIME_LIMIT) && (objective_value(model_temp) < objective_value(model_fo))
            model_fo = model_temp
        end
            
    end

    return(model_fo)
end
