# panel-method

Pressure coefficient (Cp) distribution around a cylinder using the panel method, which is a numerical technique employed to solve potential flow problems, such as the flow around a cylinder, by discretizing the surface of the object into a number of flat panels and solving for the flow parameters on each panel.

In this code, the cylinder surface is discretized into eight panels, and the coordinates of the panel vertices are calculated based on the number of panels. The panel angles, inclinations, and normal vectors are then computed. The velocity potential and the flow tangency boundary condition are utilized to derive a set of linear equations that can be solved by matrix inversion to determine the unknown panel strengths. Finally, the pressure coefficient is calculated for each panel, and the results are compared to the theoretical values of Cp for the given flow conditions.

Furthermore, the code includes some plotting commands that enable visualizing the Cp distribution and comparing it to the theoretical values.

