using PBASEX
using DelimitedFiles, HDF5

legendre_orders = collect(0:2:4)
detector_image = readdlm(joinpath(@__DIR__, "image.txt"), ' ', Float64)
println(size(detector_image))

fit = pbasex(detector_image; l=legendre_orders, basis="basis.h5")

# Julia and Python HDF5 readers may expose the matrix axes in reverse order.
h5open("out.h5", "w") do file
    println("Writing results to HDF5...")
    write(file, "E", fit.E)
    write(file, "IE", fit.IE)
    write(file, "inv_img", fit.inverted_image)
    write(file, "betas", fit.betas)
end
