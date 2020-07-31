using JuMP

function relax_and_fix(model_rf, k=2, kl=1, tlimit=600)
    t_elapsed = 0
    α = 1
    αl = -1

    nper = length(model_rf[:x][:,1])
    nprod = length(model_rf[:x][1,:])

    # relaxa as variáveis y 
    for var in model_rf[:y]
        unset_binary(var)
    end

    for i in 1:kl:nper

        # define a janela e o tempo limite de sua execução
        set_time_limit_sec(model_rf, ( (tlimit - t_elapsed) / abs((nper - α + 1) / kl) ))
        β = min(α + k -1, nper)
        if αl != -1
            α = αl
        end

        # seta como binário as variáveis y da janela
        for t in α:β
            for j in 1:nprod
                if is_fixed(model_rf[:y][t, j]) 
                    unfix(model_rf[:y][t, j])
                end
                if !is_binary(model_rf[:y][t, j]) 
                    set_binary(model_rf[:y][t, j])
                end
            end
        end

        # resolve o problema e armazena o tempo de execução
        t_elapsed += @elapsed optimize!(model_rf)

        if has_values(model_rf) && termination_status(model_rf) in (MOI.FEASIBLE_POINT, MOI.OPTIMAL, MOI.TIME_LIMIT)
            for t in α:β
                for j in 1:nprod
                    fix(model_rf[:y][t,j], value.(model_rf[:y][t,j]))
                end
            end
            α = α + kl
            αl = -1
        elseif α > 1
            αl = α - 1
        else
            error("Modelo inviável!")
        end
    
    end

    return model_rf
end