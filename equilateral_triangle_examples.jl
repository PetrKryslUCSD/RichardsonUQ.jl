using PGFPlotsX
using Statistics
using FinEtools
using Main.RichardsonExtrapolationUQ: richextrapol_uq
using CSV


# # STRI65
# relelsize = 60 ./ 2 .^(0:1:4)
# mxs = Float64[486077, 456027, 454454, 451706, 451456]
# # @show results = richextrapol_uq(mxs, relelsize)

# @show h = [r.elementsize for r in results]
# @show x = [r.estim for r in results]
# @show e = [r.estim_ad_x_2 for r in results]


# # S4
# relelsize = 60 ./ 2 .^(0:1:3)
# mxs = Float64[443426, 450060, 450892, 451455] # 451438
# results = richextrapol_uq(mxs, relelsize)

# @show h = [r.elementsize for r in results]
# @show x = [r.estim for r in results]
# @show e = [r.estim_ad_x_2 for r in results]

# STRI3
relelsize = 60 ./ 2 .^(0:1:4)
mxs = Float64[443750, 447393, 449892, 451133, 451351,  ] # ]
# relelsize = 60 ./ 2 .^(0:1:4)
# mxs = Float64[449892, 451133, 451351, 451366, 451382 ] # ]
# relelsize = 60 ./ 2 .^(0:1:2)
# mxs = Float64[451351, 451366, 451382 ] # ]
results = richextrapol_uq(mxs, relelsize)

@show h = [r.elementsize/60 for r in results]
@show x = [r.estim for r in results]
@show e = [r.estim_ad_x_2 for r in results]
@show b = [r.beta for r in results]
@show eb = [r.beta_ad_x_2 for r in results]

plots = []
@pgf p = Plot(
{
"only marks",
"error bars/y dir=both",
"error bars/y explicit",
},
Coordinates(h, x; yerror = e)
)
push!(plots, p)
@pgf p = Plot(
{
},
Coordinates([0.0, 0.275], [451389, 451389])
)
push!(plots, p)
@pgf ax = Axis(
{
"enlarge x limits=0.02",
grid = "major",
xlabel = "Relative Element Size [ND]",
ylabel = "Bending moment \$M_x\$ [N mm/mm]",
}, 
plots...
);


display(ax)



using BlackBoxOptim

# STRI3
relelsize = 60 ./ 2 .^(0:1:6)
mxs = Float64[443750, 447393, 449892, 451133, 451351, 451366, 451382 ] # ]

hcoor = [v for v in zip(relelsize, mxs)]
plots = []

@pgf p = Plot(
    {
        color = "black",
        mark = "x"
    },
    Coordinates(hcoor)
)
push!(plots, p)

# Maximum norm
result = bboptimize(x -> maximum(@. abs(x[1] - mxs - x[2] * relelsize^x[3])); SearchRange = [(0.0, 0.5e6), (0.0, 10.0), (0.1, 3.0)])
@show qex, C, beta =  best_candidate(result)

emodel(g) = qex - C * g^beta
coor = [(g, emodel(g)) for g in range(0.0, relelsize[1], length = 100)]
hcoor = [t for t in zip(relelsize, mxs)]

@pgf p = Plot(
    {
        color = "red",
    },
    Coordinates(coor)
)
push!(plots, p)

# 1-norm
result = bboptimize(x -> sum(@. abs(x[1] - mxs - x[2] * relelsize^x[3])); SearchRange = [(0.0, 0.5e6), (0.0, 10.0), (0.1, 3.0)])
@show qex, C, beta =  best_candidate(result)

emodel(g) = qex - C * g^beta
coor = [(g, emodel(g)) for g in range(0.0, relelsize[1], length = 100)]
hcoor = [t for t in zip(relelsize, mxs)]

@pgf p = Plot(
    {
        color = "green",
    },
    Coordinates(coor)
)
push!(plots, p)

# 2-norm
result = bboptimize(x -> sqrt(sum(@. abs(x[1] - mxs - x[2] * relelsize^x[3])^2)); SearchRange = [(0.0, 0.5e6), (0.0, 10.0), (0.1, 3.0)])
@show qex, C, beta =  best_candidate(result)

emodel(g) = qex - C * g^beta
coor = [(g, emodel(g)) for g in range(0.0, relelsize[1], length = 100)]
hcoor = [t for t in zip(relelsize, mxs)]

@pgf p = Plot(
    {
        color = "blue",
    },
    Coordinates(coor)
)
push!(plots, p)

@pgf ax = Axis(
    {
        xlabel = "Element size",
        ylabel = "Error"
    },
    plots...
)
display(ax)





# STRI65

println("STRI65 ========================================")


using BlackBoxOptim

relelsize_all = [
    50 ./ 2 .^(0:1:3), 
    50 ./ 2 .^(0:1:4), 
    50 ./ 2 .^(0:1:3), 
    ]
mxs_all = [
    Float64[450601, 454784, 453153, 451558,  ], 
    Float64[450601, 454784, 453153, 451558, 451439 ], 
    Float64[454784, 453153, 451558, 451439 ], 
    ]

# results = richextrapol_uq(mxs, relelsize)
let
    for data in zip(relelsize_all, mxs_all)

        println("Dataset ========================================")
        @show relelsize, mxs = data

        @show results = richextrapol_uq(mxs, relelsize)

        hcoor = [v for v in zip(relelsize, mxs)]
        plots = []

        @pgf p = Plot(
        {
        color = "black",
        mark = "x",
        "only marks",
        "thick"
        },
        Coordinates(hcoor)
        )
        push!(plots, p)

# Maximum norm
        result = bboptimize(x -> maximum(@. abs(x[1] - mxs - x[2] * relelsize^x[3]));  SearchRange = [(0.0, 0.5e6), (-15000.0, 10.0), (0.01, 3.0)]);
        @show qex, C, beta =  best_candidate(result)

        emodel(g) = qex - C * g^beta
        coor = [(g, emodel(g)) for g in range(0.0, relelsize[1], length = 100)]
        hcoor = [t for t in zip(relelsize, mxs)]

        @pgf p = Plot(
        {
        color = "red",
        style = "dashed",
        "thick"
        },
        Coordinates(coor)
        )
        push!(plots, p)

# 1-norm
        result = bboptimize(x -> sum(@. abs(x[1] - mxs - x[2] * relelsize^x[3])); SearchRange = [(0.0, 0.5e6), (-15000.0, 10.0), (0.01, 3.0)]);
        @show qex, C, beta =  best_candidate(result)

        emodel(g) = qex - C * g^beta
        coor = [(g, emodel(g)) for g in range(0.0, relelsize[1], length = 100)]
        hcoor = [t for t in zip(relelsize, mxs)]

        @pgf p = Plot(
        {
        color = "green",
        style = "dotted",
        "thick"
        },
        Coordinates(coor)
        )
        push!(plots, p)

# 2-norm
        result = bboptimize(x -> sqrt(sum(@. abs(x[1] - mxs - x[2] * relelsize^x[3])^2)); SearchRange = [(0.0, 0.5e6), (-15000.0, 10.0), (0.01, 3.0)]);
        @show qex, C, beta =  best_candidate(result)

        emodel(g) = qex - C * g^beta
        coor = [(g, emodel(g)) for g in range(0.0, relelsize[1], length = 100)]
        hcoor = [t for t in zip(relelsize, mxs)]

        @pgf p = Plot(
        {
        color = "blue",
        style = "dashdotted",
        "thick"
        },
        Coordinates(coor)
        )
        push!(plots, p)


        @pgf ax = Axis(
        {
        xlabel = "Element size",
        ylabel = "Error",
        grid = "major",
        },
        plots...
        )
        display(ax)

    end
end

nothing