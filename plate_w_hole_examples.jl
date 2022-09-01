module plate_w_hole_examples
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsDeforLinear.AlgoDeforLinearModule
using FinEtools.MeshExportModule
using LinearAlgebra: norm
using SparseArrays: cholesky
using Statistics: mean
using AlgebraicMultigrid, LinearOperators, Krylov
using LinearAlgebra, SparseArrays, MKLSparse, Metis

function amg(A::AbstractMatrix{T}; kwargs...) where T
  n, m = size(A)
  ml = ruge_stuben(A; kwargs...)
  # ml = smoothed_aggregation(A; kwargs...)
  P = aspreconditioner(ml)
  return LinearOperator(T, n, n, true, true, (y, v) -> ldiv!(y, P, v))
end

atol =  0.0
rtol = 1.0e-5

function plate_w_hole_H20_stress(alpha = 0.7, nref = 6)
    E = 2.4*phun("MEGA*PA");# 210e3 MPa
    nu = 0.49995;
    Re = 0.3*phun("M"); # outer radius
    Ri= 0.1*phun("M"); # hole radius
    H = 0.1*phun("M") # thickness of the plate
    nRadial, nCircumferential, nThickness = 6, 8, 1;
    sigma0=1*phun("MEGA*PA");
    CSVFile = "plate_w_hole_H20_stress.CSV"

    function sigmaxx(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r=norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r=norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end

    # convergencestudy = FDataDict[]
    sigxsb = FFlt[]
    sigysa = FFlt[]
    numelements = Int[]
    numnodes = Int[]
    relelsize = FFlt[]
    previous_nmult = -2
    for ref in 0:1:nref
        println("ref = $(ref)")
        nmult = Int(round(1/alpha^ref))
        nmult <= previous_nmult+2  && (nmult = previous_nmult+2)
        @show previous_nmult = nmult
        Thickness = H
        tolerance = Thickness/nmult/1000.; # Geometrical tolerance
        
        fens,fes = H20block(1.0, pi/2, Thickness, nmult*nRadial, nmult*nCircumferential, nmult*nThickness)
        @show count(fens)

        bdryfes = meshboundary(fes);
        icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);

        for i=1:count(fens)
            t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
            fens.xyz[i,:] = [(t*Re+(1-t)*Ri)*cos(a), (t*Re+(1-t)*Ri)*sin(a), z];
        end

        geom = NodalField(fens.xyz)
        u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

        l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
        setebc!(u,l1,true, 2, 0.0)
        l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
        setebc!(u,l1,true, 1, 0.0)
        # PLANE-STRESS constraint: assume the plane z=0 is the plane of symmetry of the plate
        l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
        setebc!(u,l1,true, 3, 0.0)
        # If, in addition to the above, this was enabled, the PLANE-STRAIN
        # constraint would be enforced.
        l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
        setebc!(u,l1,true, 3, 0.0)

        applyebc!(u)
        numberdofs!(u)
        el1femm =  FEMMBase(IntegDomain(subset(bdryfes,icl), GaussRule(2, 3)))
        function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
            local r = sqrt(XYZ[1]^2 + XYZ[2]^2)
            nx = XYZ[1]/r; ny = XYZ[2]/r
            forceout[1] = sigmarr(XYZ) * nx - sigmart(XYZ) * ny
            forceout[2] = sigmarr(XYZ) * ny + sigmart(XYZ) * nx
            forceout[3] = 0.0
            return forceout
        end
        fi = ForceIntensity(FFlt, 3, pfun);
        F2 = distribloads(el1femm, geom, u, fi, 2);

        MR = DeforModelRed3D

        material = MatDeforElastIso(MR, E, nu)

        femm = FEMMDeforLinear(MR, IntegDomain(fes, GaussRule(3, 2)), material)

        # The geometry field now needs to be associated with the FEMM
        femm = associategeometry!(femm, geom)

        K = stiffness(femm, geom, u)
        # K = cholesky(K)
        # U = K\(F2)
        P = amg(K)
        U, stats = minres_qlp(K, F2, M=P, atol=atol, rtol=rtol, verbose=0)
        scattersysvec!(u, U[:])

        nlA = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
        nlB = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, Thickness], inflate=tolerance);
        
        sigx = fieldfromintegpoints(femm, geom, u, :Cauchy, 1)
        sigy = fieldfromintegpoints(femm, geom, u, :Cauchy, 2)
        sigyA = mean(sigy.values[nlA,1], dims = 1)[1]
        sigyAtrue = sigmatt([Ri, 0.0, 0.0])
        println("sig_y@A =$(sigyA/phun("MPa")) vs $(sigyAtrue/phun("MPa")) [MPa]")
        sigxB = mean(sigx.values[nlB,1], dims = 1)[1]
        sigxBtrue = sigmatt([0.0, Ri, 0.0])
        println("sig_x@B =$(sigxB/phun("MPa")) vs $(sigxBtrue/phun("MPa")) [MPa]")
        push!(numnodes, count(fens))
        push!(numelements, count(fes))
        push!(relelsize, alpha^ref)
        push!(sigxsb, sigxB)
        push!(sigysa, sigyA)

        File =  "plate_w_hole_H20_stress-sigx -$ref.vtk"
        vtkexportmesh(File, fes.conn, geom.values,
        FinEtools.MeshExportModule.VTK.H20; vectors=[("u", u.values)],        scalars=[("sigmax", sigx.values/phun("MEGA*PA")), ("sigmay", sigy.values/phun("MEGA*PA"))])
        
        savecsv(CSVFile, numelements=vec(numelements), numnodes=vec(numnodes), 
            relelsize=vec(relelsize), sigxsb=vec(sigxsb), sigysa=vec(sigysa))
        
    end # for ref in


    return numnodes, numelements, relelsize, sigxsb, sigysa
end # plate_w_hole_H20_stress


function plate_w_hole_MSH8_stress(nu = 0.4999, alpha = 0.5, nref = 5)
    E = 2.4*phun("MEGA*PA");# 210e3 MPa
    Re = 0.3*phun("M"); # outer radius
    Ri= 0.1*phun("M"); # hole radius
    H = 0.1*phun("M") # thickness of the plate
    nRadial, nCircumferential, nThickness = 4, 5, 1;
    sigma0=1*phun("MEGA*PA");
    CSVFile = "plate_w_hole_MSH8_stress.CSV"

    function sigmaxx(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0*(1-Ri^2/r^2*(3/2*cos(2*th)+cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmayy(x)
        local r=norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*cos(2*th)-cos(4*th))+3/2*Ri^4/r^4*cos(4*th));
    end
    function sigmaxy(x)
        local r=norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0*(Ri^2/r^2*(1/2*sin(2*th)+sin(4*th))-3/2*Ri^4/r^4*sin(4*th));
    end
    function sigmarr(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0/2*(1-Ri^2/r^2) + sigma0/2*(1-4*Ri^2/r^2+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmatt(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return sigma0/2*(1+Ri^2/r^2) - sigma0/2*(1+3*Ri^4/r^4)*cos(2*th)
    end
    function sigmart(x)
        local r = norm(x[1:2]);
        local th = atan(x[2],x[1]);
        return -sigma0/2*(1+2*Ri^2/r^2-3*Ri^4/r^4)*sin(2*th)
    end

    # convergencestudy = FDataDict[]
    sigxsb = FFlt[]
    sigysa = FFlt[]
    numelements = Int[]
    numnodes = Int[]
    relelsize = FFlt[]
    previous_nmult = -2
    for ref in 0:1:nref
        println("ref = $(ref)")
        nmult = Int(round(1/alpha^ref))
        # nmult <= previous_nmult+2  && (nmult = previous_nmult+2)
        # @show previous_nmult = nmult
        Thickness = H
        tolerance = Thickness/nmult/1000.; # Geometrical tolerance
        
        fens,fes = H8block(1.0, pi/2, Thickness, nmult*nRadial, nmult*nCircumferential, nmult*nThickness)
        @show count(fens)

        bdryfes = meshboundary(fes);
        icl = selectelem(fens, bdryfes, box=[1.0, 1.0, 0.0, pi/2, 0.0, Thickness], inflate=tolerance);

        for i=1:count(fens)
            t=fens.xyz[i,1]; a=fens.xyz[i,2]; z=fens.xyz[i,3]
            fens.xyz[i,:] = [(t*Re+(1-t)*Ri)*cos(a), (t*Re+(1-t)*Ri)*sin(a), z];
        end

        geom = NodalField(fens.xyz)
        u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

        l1 =selectnode(fens; box=[0.0, Inf, 0.0, 0.0, 0.0, Thickness], inflate = tolerance)
        setebc!(u,l1,true, 2, 0.0)
        l1 =selectnode(fens; box=[0.0, 0.0, 0.0, Inf, 0.0, Thickness], inflate = tolerance)
        setebc!(u,l1,true, 1, 0.0)
        # PLANE-STRESS constraint: assume the plane z=0 is the plane of symmetry of the plate
        l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, 0.0, 0.0], inflate = tolerance)
        setebc!(u,l1,true, 3, 0.0)
        # If, in addition to the above, this was enabled, the PLANE-STRAIN
        # constraint would be enforced.
        l1 =selectnode(fens; box=[0.0, Inf, 0.0, Inf, Thickness, Thickness], inflate = tolerance)
        setebc!(u,l1,true, 3, 0.0)

        applyebc!(u)
        numberdofs!(u)
        el1femm =  FEMMBase(IntegDomain(subset(bdryfes,icl), GaussRule(2, 3)))
        function pfun(forceout::FVec{T}, XYZ::FFltMat, tangents::FFltMat, fe_label::FInt) where {T}
            local r = sqrt(XYZ[1]^2 + XYZ[2]^2)
            nx = XYZ[1]/r; ny = XYZ[2]/r
            forceout[1] = sigmarr(XYZ) * nx - sigmart(XYZ) * ny
            forceout[2] = sigmarr(XYZ) * ny + sigmart(XYZ) * nx
            forceout[3] = 0.0
            return forceout
        end
        fi = ForceIntensity(FFlt, 3, pfun);
        F2 = distribloads(el1femm, geom, u, fi, 2);

        MR = DeforModelRed3D

        material = MatDeforElastIso(MR, E, nu)

        femm = FEMMDeforLinearMSH8(MR, IntegDomain(fes, GaussRule(3, 2)), material)

        # The geometry field now needs to be associated with the FEMM
        femm = associategeometry!(femm, geom)

        K = stiffness(femm, geom, u)
        # K = cholesky(K)
        # U = K\(F2)
        P = amg(K)
        U, stats = minres_qlp(K, F2, M=P, atol=atol, rtol=rtol, verbose=0)
        scattersysvec!(u, U[:])

        nlA = selectnode(fens, box=[Ri, Ri, 0.0, 0.0, 0.0, Thickness], inflate=tolerance);
        nlB = selectnode(fens, box=[0.0, 0.0, Ri, Ri, 0.0, Thickness], inflate=tolerance);
        
        sigx = fieldfromintegpoints(femm, geom, u, :Cauchy, 1)
        sigy = fieldfromintegpoints(femm, geom, u, :Cauchy, 2)
        sigyA = mean(sigy.values[nlA,1], dims = 1)[1]
        sigyAtrue = sigmatt([Ri, 0.0, 0.0])
        println("sig_y@A =$(sigyA/phun("MPa")) vs $(sigyAtrue/phun("MPa")) [MPa]")
        sigxB = mean(sigx.values[nlB,1], dims = 1)[1]
        sigxBtrue = sigmatt([0.0, Ri, 0.0])
        println("sig_x@B =$(sigxB/phun("MPa")) vs $(sigxBtrue/phun("MPa")) [MPa]")
        push!(numnodes, count(fens))
        push!(numelements, count(fes))
        push!(relelsize, alpha^ref)
        push!(sigxsb, sigxB)
        push!(sigysa, sigyA)

        File =  "plate_w_hole_MSH8_stress-sigx -$ref.vtk"
        vtkexportmesh(File, fes.conn, geom.values,
        FinEtools.MeshExportModule.VTK.H8; vectors=[("u", u.values)],        scalars=[("sigmax", sigx.values/phun("MEGA*PA")), ("sigmay", sigy.values/phun("MEGA*PA"))])
        
        savecsv(CSVFile, numelements=vec(numelements), numnodes=vec(numnodes), 
            relelsize=vec(relelsize), sigxsb=vec(sigxsb), sigysa=vec(sigysa))
        
    end # for ref in


    return numnodes, numelements, relelsize, sigxsb, sigysa
end # plate_w_hole_H20_stress



# numnodes, numelements, relelsize, sigxsb, sigysa = plate_w_hole_MSH8_stress()

using PGFPlotsX
using Statistics
using FinEtools
using FinEtools.AlgoBaseModule: richextrapol
using Main.RichardsonExtrapolationUQ: richextrapol_uq
using CSV

f = CSV.File("plate_w_hole_MSH8_stress.CSV")     
relelsize = [l.relelsize for l in f]
sigxsb = [l.sigxsb for l in f]
@show results = richextrapol_uq(sigxsb, relelsize)

h = [r.elementsize for r in results]
@show x = [r.estim ./ phun("MPa") for r in results]
e = [r.estim_ad_x_2 ./ phun("MPa") for r in results]
plots = []
@pgf p = Plot(
{
"only marks",
"error bars/y dir=both",
"error bars/y explicit",
"tick style=thick"
},
Coordinates(h, x; yerror = e)
)
push!(plots, p)
@pgf p = Plot(
{
"thick"
},
Coordinates([0.0, 0.25], [3.0, 3.0])
)
push!(plots, p)
@pgf ax = Axis(
{
height = "7cm",
width = "7cm",
"enlarge x limits=0.02",
ymin = 1.75,
ymax = 3.5,
grid = "major",
xlabel = "Relative element size [ND]",
ylabel = "Stress \$\\sigma_x\$ [MPa]",
scaled_ticks = false
}, 
plots...
);

display(ax)


using PGFPlotsX
using Statistics
using FinEtools
using FinEtools.AlgoBaseModule: richextrapol
using Main.RichardsonExtrapolationUQ: gci_uq
using CSV

f = CSV.File("plate_w_hole_MSH8_stress.CSV")     
relelsize = [l.relelsize for l in f]
sigxsb = [l.sigxsb for l in f]
@show results = gci_uq(sigxsb, relelsize)

h = [r.elementsize for r in results]
x = [r.w_1 ./ phun("MPa") for r in results]
e = [r.gci ./ phun("MPa") * r.w_1 for r in results]

@show [x e]

b = [r.beta for r in results]

plots = []
@pgf p = Plot(
{
"only marks",
"error bars/y dir=both",
"error bars/y explicit",
"tick style=thick"
},
Coordinates(h, x; yerror = e)
)
push!(plots, p)
@pgf p = Plot(
{
"thick"
},
Coordinates([0.0, 0.25], [3.0, 3.0])
)
push!(plots, p)
@pgf ax = Axis(
{
height = "7cm",
width = "7cm",
"enlarge x limits=0.02",
ymin = 1.75,
ymax = 3.5,
grid = "major",
xlabel = "Relative element size [ND]",
ylabel = "Stress \$\\sigma_x\$ [MPa]",
scaled_ticks = false
}, 
plots...
);

display(ax)



using PGFPlotsX
using Statistics
using FinEtools
using FinEtools.AlgoBaseModule: richextrapol
using Main.RichardsonExtrapolationUQ: bsr_uq
using CSV

f = CSV.File("plate_w_hole_MSH8_stress.CSV")     
relelsize = [l.relelsize for l in f]
sigxsb = [l.sigxsb for l in f]
@show results = bsr_uq(sigxsb, relelsize)

h = [r.elementsize for r in results]
x = [r.q_m ./ phun("MPa") for r in results]
e = [r.epshat_m * r.q_m ./ phun("MPa") for r in results]

@show [x e]

b = [r.chat_m for r in results]

plots = []
@pgf p = Plot(
{
"only marks",
"error bars/y dir=both",
"error bars/y explicit",
"tick style=thick"
},
Coordinates(h, x; yerror = e)
)
push!(plots, p)
@pgf p = Plot(
{
"thick"
},
Coordinates([0.0, 0.25], [3.0, 3.0])
)
push!(plots, p)
@pgf ax = Axis(
{
height = "7cm",
width = "7cm",
"enlarge x limits=0.02",
ymin = 1.75,
ymax = 3.5,
grid = "major",
xlabel = "Relative element size [ND]",
ylabel = "Stress \$\\sigma_x\$ [MPa]",
scaled_ticks = false
}, 
plots...
);

display(ax)

end # module plate_w_hole_examples

