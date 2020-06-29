using JuMP

function relax_and_fix(model, k=2, kl=1, tlimit=600)
    model_rf = model
    t_elapsed = 0
    α = 1
    nper = length(model[:x][:,1])
    nprod = length(model[:x][1,:])

    # relaxa as variáveis y
    for var in model_rf[:y]
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

        # define a janela e o tempo limite de sua execução
        set_time_limit_sec(model_rf, ( (tlimit - t_elapsed) / abs((nper - α + 1) / kl) ))
        β = min(α + k -1, nper)

        # seta como binário as variáveis y da janela
        for t in α:β
            for j in 1:nprod
                #set_integer(model[:x][t, j])
                set_binary(model_rf[:y][t, j])
                #set_integer(model[:s][t, j])
            end
        end

        # resolve o problema e armazena o tempo de execução
        t_elapsed += @elapsed optimize!(model_rf)

        if has_values(model_rf) && termination_status(model_rf) in (MOI.FEASIBLE_POINT, MOI.OPTIMAL, MOI.TIME_LIMIT)
            for t in α:(α + kl -1)
                for j in 1:nprod
                    #fix(model[:x][t,j], value.(model[:x][t,j]))
                    fix(model_rf[:y][t,j], value.(model_rf[:y][t,j]))
                    #fix(model[:s][t,j], value.(model[:s][t,j]))
                end
            end
            α = α + kl
            β = min(α + k -1, nper)
        else
            # TODO: trabalhar a exceção de quando a solução for inviável
            error("Status do modelo deu ruim")
        end
    
    end

    return model_rf
end