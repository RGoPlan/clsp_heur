module Instances

using DelimitedFiles
using JuMP, Gurobi

export TTMcCinstance, read_instance_TTMcC, model_TTMcC

# cria um ambiente com o Gurobi solver
print("Gurobi solver: ")
const GUROBI_ENV = Gurobi.Env()

# cria o tipo TTMcCinstance para formatar a instância
mutable struct TTMcCinstance
    nprod::Int16    
    nper::Int16
    p::Float16              # custo unitário de produção do item j
    c::Int32                # capacidade de produção + estoque
    a::Array{Float16}       # insumos consumidos para produzir uma unidade do item j
    b::Array{Float16}       # insumos consumidos para setup de uma unidade do item j
    h::Array{Float16}       # custo unitário de estocagem do item j
    q::Array{Float16}       # custo unitário de setup do item j
    d::Array{Float16}       # demanda determinística do item j no período t
end

# ler o arquivo de texto e retorna um tipo TTMcCinstance
function read_instance_TTMcC(file_path::String)
    f = readdlm(file_path)

    nprod = f[1,1]
    nper = f[1,2]

    inst = TTMcCinstance(
        f[1,1],
        f[1,2],
        f[2,1],
        f[3,1],
        f[4:4+nprod-1, 1],
        f[4:4+nprod-1, 2],
        f[4:4+nprod-1, 3],
        f[4:4+nprod-1, 4],
        f[4+nprod:4+nprod+nper-1, 1:nprod]
    )

    return inst
end

# recebe um TTMcCinstance e modela o problema
function model_TTMcC(inst::TTMcCinstance)
    model = Model(() -> Gurobi.Optimizer(GUROBI_ENV))

    @variables(model, begin
        x[t=1:inst.nper, j=1:inst.nprod], Int
        y[t=1:inst.nper, j=1:inst.nprod], Bin
        s[t=1:inst.nper, j=1:inst.nprod], Int
    end)

    @objective(
        model, Min,
        sum(inst.p*x[t,j] + inst.q[j]*y[t,j] + inst.h[j]*s[t, j] for t=1:inst.nper, j=1:inst.nprod)
    )

    @constraints(model, begin
        capac_constr[t=1:inst.nper], sum(inst.a[j]*x[t,j] + inst.b[j]*y[t,j] for j=1:inst.nprod) <= inst.c
        stock1[j=1:inst.nprod], x[1, j] - s[1, j] == inst.d[1, j]
        stock[t=2:inst.nper, j=1:inst.nprod], inst.d[t, j] + s[t, j] == x[t, j] + s[t-1,j]
        setup[t=1:inst.nper, j=1:inst.nprod], x[t,j] <= y[t,j]*max((inst.c - inst.b[j]) / inst.a[j], sum(inst.d[t,j] for t=t:inst.nper))
        nnegatx[t=1:inst.nper, j=1:inst.nprod], x[t,j] >= 0
        nnegats[t=1:inst.nper, j=1:inst.nprod], s[t,j] >= 0
    end)

    return model
end

end