# MeepMedium
Python code to calculate permittivity from Meep medium

Meep is an electromagnetic simulation package based on FDTD method. We have to define the permittivity and permeability of all kinds of materials in the simulations. But the definition process is quite indirect. We have to define a series of Lorentz or Drude resonators. Then the permittivity and permeability can be calculated through a complex equation. So here I wrote several functions to calculte the permittivty from Meep medium class.
