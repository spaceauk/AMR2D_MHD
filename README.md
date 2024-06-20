# 2D MHD Riemann solver with AMR 
- AMR is based on a quadtree data structure. The coarsening and refining thresholds used are 0.05 and 1.0, respectively.
- Contains constrained transport to ensure solenoidal magnetic field.
- Generalized boundary conditions for higher order schemes like 5th-order WENO method.

## Sample plots
Following plots are started with 4 root blocks each with grid size of 16 by 16. 
1) 1D Brio-Wu shock tube at t=0.2. 4 levels of grid refinement are used.

   a. Density:
   ![BWx022t_0 200369_16x16_MC_RUSA_density](https://github.com/spaceauk/AMR2D_MHD/assets/64028216/c5b751c8-485d-4b9e-83a0-38ea50a4f447)


3) 2D MHD rotor at t=0.295. 5 levels of grid refinement are used.
   
   a. Density:
   ![rotor045t_0 295062_16x16_lv5_MC_RUSA_density](https://github.com/spaceauk/AMR2D_MHD/assets/64028216/231fdce4-20ae-4040-b503-fe74af8072dc)
   b. Pressure:
file:///home/spaceghoul/Desktop/Summer2024/test_AMR/MHD_MUSCL_n_WENO_genBC_n_AMR/plots/good_img_AMR/rotor045t_0.295062_16x16_lv5_MC_RUSA_pressure.png
   c. Magnetic field along x:  file:///home/spaceghoul/Desktop/Summer2024/test_AMR/MHD_MUSCL_n_WENO_genBC_n_AMR/plots/good_img_AMR/rotor045t_0.295062_16x16_lv5_MC_RUSA_Bx.png



