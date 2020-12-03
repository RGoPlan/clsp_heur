using JuMP

function relax_and_fix(model_rf, k=2, kl=1, tlimit=120)
    t_elapsed = 0
    α = 1
    αl = -1

    nper = length(model_rf[:x][:,1])
    nprod = length(model_rf[:x][1,:])

    β = min(α + k -1, nper)

    # relaxa as variáveis y 
    for var in model_rf[:y]
        unset_binary(var)
    end

    n_inter = 0
    n_bt = 0
    while β < nper
        # define a janela e o tempo limite de sua execução
        set_time_limit_sec(model_rf, ( (tlimit - t_elapsed) / abs((nper - α + 1) / kl) ))
        #println("tempo setado")
        
        if αl != -1
            α = αl
        else
            β = min(α + k -1, nper)
        end

        #println("[$(α), $(β)]")

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
            # println("Estou voltando...")
            n_bt += 1
            αl = α - 1
            if β == 20
                β -= 1
            end
        else
            break
        end
        n_inter += 1
    
    end

    return model_rf, n_inter, n_bt, t_elapsed
end