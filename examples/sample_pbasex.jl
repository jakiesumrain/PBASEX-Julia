include("../src/pbasex_core.jl")
using .pBASEX
using DelimitedFiles, HDF5

# Define the basis function parameters
l_values = collect(0:2:4)  # Legendre polynomial orders


img = readdlm("image.txt", ' ', Float64)
println(size(img))

out = pbasex(img, l = l_values, basis = "basis.h5")

# Be aware that when you read a h5 file saved in python with julia, the matrices may be transposed. The same is true vise visa.
h5open("out.h5", "w") do file
    println("Writing results to HDF5...") # Changed print statement for clarity
    write(file, "E", out.E)
    write(file, "IE", out.IE)
    write(file, "inv_img", out.inverted_image)
    write(file, "betas", out.betas)
end

#println(minimum(basis_file.S))
#println(maximum(basis_file.S))

# add a histogram of the singular values
#histogram(basis_file.S)
#savefig("singular_values_histogram_1.png")

# Save singular values to text file
#open("singular_values_1.txt", "w") do io
    #writedlm(io, basis_file.S)
#end

# add a plot of the singular values
#plot(basis_file.S)
#savefig("singular_values_plot_1.png")


