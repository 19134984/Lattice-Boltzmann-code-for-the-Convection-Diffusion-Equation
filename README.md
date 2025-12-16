# Lattice-Boltzmann-code-for-the-Convection-Diffusion-Equation
Solving linear and nonlinear convection-diffusion equations using the Lattice Boltzmann Method with MRT D2Q5 or D2Q9 models.
Gaussianhill文件是线性对流扩散方程，无源项，各向同性的高斯山扩散，其中由于宏观方程是有量纲的，首先通过单位转换到格子单位下，除了被动标量Φ以外，其余均需要单位转换，然后计算时间选择了物理时间的1s，由此计算需要多少格子时间，主要参考文献为Huang, Rongzong, and Huiying Wu. 2014. “A Modified Multiple-Relaxation-Time Lattice Boltzmann Model for Convection–Diffusion Equation.” Journal of Computational Physics 274:50–63. https://doi.org/10.1016/j.jcp.2014.05.041.


Gaussianhillanisotropic文件是线性对流扩散方程，无源项，各向异性的高斯山扩散，首先也得单位换算，但由于扩散系数是一个矩阵，于是得选取一个物理扩散张量的最大特征值 λ_max_phys来代表扩散强度，并借此计算1单位格子时间代表多少物理时间，高斯山这两个算例均是考虑了周期边界（按照其他文献所说，扩散在边界很小，可以视为周期边界），主要参考的还是文献Huang, Rongzong, and Huiying Wu. 2014. “A Modified Multiple-Relaxation-Time Lattice Boltzmann Model for Convection–Diffusion Equation.” Journal of Computational Physics 274:50–63. https://doi.org/10.1016/j.jcp.2014.05.041.


linearADEwithasource文件是线性对流扩散方程，有源项，各向同性的基准算例，主要来自于文献的算例3.1：Guo, Xiuya, Baochang Shi, and Zhenhua Chai. 2018. “General Propagation Lattice Boltzmann Model for Nonlinear Advection-Diffusion Equations.” Physical Review E 97 (4). https://doi.org/10.1103/PhysRevE.97.043310.
以及文献的例1：Shang, Jinlong, Zhenhua Chai, Huili Wang, and Baochang Shi. 2020. “Discrete Unified Gas Kinetic Scheme for Nonlinear Convection-Diffusion Equations.” Physical Review E 101 (2). https://doi.org/10.1103/PhysRevE.101.023306.


nonlinearheatconduction文件是非线性扩散方程（无对流），有源项，各向同性的基准算例，主要来自于献的算例3.2：Guo, Xiuya, Baochang Shi, and Zhenhua Chai. 2018. “General Propagation Lattice Boltzmann Model for Nonlinear Advection-Diffusion Equations.” Physical Review E 97 (4). https://doi.org/10.1103/PhysRevE.97.043310.
以及文献的例3：Shang, Jinlong, Zhenhua Chai, Huili Wang, and Baochang Shi. 2020. “Discrete Unified Gas Kinetic Scheme for Nonlinear Convection-Diffusion Equations.” Physical Review E 101 (2). https://doi.org/10.1103/PhysRevE.101.023306.
