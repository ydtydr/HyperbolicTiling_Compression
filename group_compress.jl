include("./Hyperbolic.jl")
using Main.Hyperbolic
using StaticArrays
using LinearAlgebra
using DelimitedFiles
using ArgParse

argparse = ArgParseSettings()
@add_arg_table argparse begin
    "--dataset"
        help = "which dataset to compress: gr-qc, wordnet or bio-yeast"
        arg_type = String
        required = true
    "--compr_dir"
        help = "file directory to store the compressed L-tiling model, txt format"
        arg_type = String
        required = true
    "--model"
        help = "what models to be compressed, poincare or lorentz"
        arg_type = String
        default = "poincare"
end
args = parse_args(argparse)

if args.dataset == "gr-qc"
    preci = 500
    node_num = 4158
    embedding_file = "./dataset/grqc_500.emb"
elseif args.dataset == "wordnet"
    preci = 8200
    node_num = 74374
    embedding_file = "./dataset/wordnet2h-comb.emb"
    if args.model == "lorentz"
        embedding_file = "wordnet_lorendim3.txt"
    end
elseif args.dataset == "bio-yeast"
    preci = 4000
    node_num = 1458
    embedding_file = "./dataset/bio-yeast.r2.emb"
end

setprecision(preci)
###################
#To load a poincare embedding of Gr-Qc dataset
poincare_pts = zeros(BigFloat, node_num+1, 2);
totaltime, totallines = open(embedding_file) do f
    linecounter = 0
    timetaken = @elapsed for li in eachline(f)
        linecounter += 1
        firs = findfirst(isequal(','), li)
        secon = findnext(isequal(','), li,firs+1)
        thir = findnext(isequal(','), li,secon+1)
        poincare_pts[linecounter,1] = BigFloat(li[firs+1:secon-1])
        if show(thi)==nothing
            poincare_pts[linecounter,2] = BigFloat(li[secon+1:end])
        else
            poincare_pts[linecounter,2] = BigFloat(li[secon+1:thir-1])
        end
    end
    println(linecounter)
    (timetaken, linecounter)
end
poincare_pts = poincare_pts[2:end,:];
#############################
##or you can directly load a lorentz embedding like this
# node_num = 4158
# lorentz_pts = zeros(BigFloat, node_num, 3);
# totaltime, totallines = open(embedding_file) do f
#     linecounter = 0
#     timetaken = @elapsed for li in eachline(f)
#         linecounter += 1
#         fir = findfirst(isequal(','), li)
#         sec = findnext(isequal(','), li,fir+1)
#         lorentz_pts[linecounter,1] = BigFloat(li[1:fir-1])
#         lorentz_pts[linecounter,2] = BigFloat(li[fir+1:sec-1])
#         lorentz_pts[linecounter,3] = BigFloat(li[sec+1:end])
#     end
#     (timetaken, linecounter)
# end

let
##To get L-tiling model points
#poincare_pts or lorentz_pts should be bigfloat matrices
out = open(compr_dir,"w")# file to store VBW encoding the points in the fundamental domain
erro = 0
@time begin # calculate time
    for num1=1:node_num
        # if load poincare_pts
        hyp_poincare_pt = Main.Hyperbolic.HyPPoincare(SVector{2}(poincare_pts[num1,:]))
        ###transform poincare embedding to lorentz embedding
        hyp_lorentz_pt = Main.Hyperbolic.HyPLorentz(hyp_poincare_pt)
        # if load lorentz_pts
#         hyp_lorentz_pt = Main.Hyperbolic.HyPLorentz(SVector{3}(lorentz_pts[num1,:]))
        setprecision(preci)
        fpt_vbw = Main.Hyperbolic.Gen_VBW(hyp_lorentz_pt, BigInt);
        #map from L-tiling_vbw to lorenze model so as to calculate error, this step may cost some time
        hyp_lorentz_pt_re = Main.Hyperbolic.HyPLorentz(fpt_vbw)
        erro += Main.Hyperbolic.dist(hyp_lorentz_pt,hyp_lorentz_pt_re)^2;
        write(string(fpt_vbw.c[1]),",",string(fpt_vbw.c[2]),",",string(fpt_vbw.c[3]), fpt_vbw.vbw,"\n");
    end
    close(out)
    avg_error = Float64(sqrt(erro/node_num))
    println("Average mean squared hyperbolic distance is ", avg_error)
end
end
println("Model Compressed Successfully")