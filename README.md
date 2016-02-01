# ns3D
Three-Dimensional Mixed-Element Unstructured Finite-Volume Solver for Compressible Euler/Navier-Stokes Equations

#Quick Intro
1- Open a terminal and go to ./src/ folder and type "make"

2- execute the generated binary "ns3d". The code starts to run a compressible viscous Blasius boundary layer over a flat plate. It finally converges to the solution below

![alt tag](https://raw.github.com/arrgasm/ns3D/master/imgs/uCont.png)
![alt tag](https://raw.github.com/arrgasm/ns3D/master/imgs/uVec.png)

The solution is compared to the Blasius exact solution and a solution obtained by commercial solver Ansys-Fluent on the same grid. The result is presented below.

![alt tag](https://raw.github.com/arrgasm/ns3D/master/imgs/Compare.png)

This code also captures crisp shocks. See the following supersonic ramp! At this point, you can go ahead and read main.c file and modify/extend it for your own work. Good luck!

![alt tag](https://raw.github.com/arrgasm/ns3D/master/imgs/composite3D.png)


