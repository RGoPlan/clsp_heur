using JuMP

function fix_and_optimize(model_fo, k=2, kl=1, tlimit=600)
    t_elapsed = 0
    obj_iter = objective_value(model_fo)
    last_obj_iter = Inf

    nper = length(model_fo[:x][:,1])
    nprod = length(model_fo[:x][1,:])

    set_time_limit_sec(model_fo, ( (tlimit - t_elapsed) / (nper / kl) ))

    while (time_limit_sec(model_fo)) > 0 && (obj_iter < last_obj_iter)
        last_obj_iter = obj_iter
        α = 1

        for i in range(1, stop=nper, step=kl)
            β = min(α + k -1, nper)

            # retira os valores fixos das variáveis y da janela
            for t in α:β
                for j in 1:nprod
                    unfix(model_fo[:y][t,j])
                end
            end

            # resolve o problema e armazena o tempo de execução
            t_elapsed = @elapsed optimize!(model_fo)
            set_time_limit_sec(model_fo, ( (tlimit - t_elapsed) / (abs(nper - α + 1) / kl) ))

            for t in α:β
                for j in 1:nprod
                    fix(model_fo[:y][t,j], value.(model_fo[:y][t,j]))
                end
            end

            α = α + kl

        end

        obj_iter = objective_value(model_fo)
    end

    return(model_fo)
end
