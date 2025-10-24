# setup paths -
const _ROOT = @__DIR__;
const _PATH_TO_DATA = joinpath(_ROOT, "data");
const _PATH_TO_SRC = joinpath(_ROOT, "src");

# if we are missing any packages, install them -
using Pkg;
if (isfile(joinpath(_ROOT, "Manifest.toml")) == false) # have manifest file, we are good. Otherwise, we need to instantiate the environment
    Pkg.add(path="https://github.com/varnerlab/VLDataScienceMachineLearningPackage.jl.git")
    Pkg.activate("."); Pkg.resolve(); Pkg.instantiate(); Pkg.update();
end

# load external packages -
using VLDataScienceMachineLearningPackage
using LinearAlgebra
using BenchmarkTools
using Statistics
using Test
using Images
using TestImages
using ImageMagick
using ImageIO
using DelimitedFiles
using Plots
using DataFrames
using Random
using Distributions
using PrettyTables
using JuMP
using MadNLP
using MathOptInterface

# setup random number generator -
Random.seed!(1234); # seed the random number generator for reproducibility

# include our source files -
include(joinpath(_PATH_TO_SRC, "Types.jl"));
include(joinpath(_PATH_TO_SRC, "Factory.jl"));
include(joinpath(_PATH_TO_SRC, "Compute.jl"));
