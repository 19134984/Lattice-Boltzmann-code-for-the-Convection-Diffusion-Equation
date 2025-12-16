# Lattice-Boltzmann-code-for-the-Convection-Diffusion-Equation
Solving linear and nonlinear convection-diffusion equations using the Lattice Boltzmann Method with MRT D2Q5 or D2Q9 models.
gaosishan文件是线性对流扩散方程，无源项，各向同性的高斯山扩散，其中由于宏观方程是有量纲的，首先通过单位转换到格子单位下，除了被动标量Φ以外，其余均需要单位转换，然后计算时间选择了物理时间的1s，由此计算需要多少格子时间，主要参考文献为Huang, Rongzong, and Huiying Wu. 2014. “A Modified Multiple-Relaxation-Time Lattice Boltzmann Model for Convection–Diffusion Equation.” Journal of Computational Physics 274:50–63. https://doi.org/10.1016/j.jcp.2014.05.041.
