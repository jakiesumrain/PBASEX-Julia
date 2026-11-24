using pBASEX

detector_coordinates = collect(0:511)
legendre_orders = collect(0:2:4)

# Save the projected basis and central-slice basis for reuse.
basis_data = BASISData{Float64}("./basis.h5"; x=detector_coordinates, l=legendre_orders)
