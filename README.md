# 2D MHD Riemann solver with AMR 
- AMR is based on a quadtree data structure. The coarsening and refining thresholds used are 0.05 and 1.0, respectively.
- Contains constrained transport to ensure solenoidal magnetic field.
- Generalized boundary conditions for higher order schemes like 5th-order WENO method.

## Sample plots
Following plots are started with 4 root blocks each with grid size of 16 by 16. 
1) 1D Brio-Wu shock tube at t=0.2. 4 levels of grid refinement are used.

   a. Density:
   ![BWx022t_0 200369_16x16_MC_RUSA_density](https://github.com/spaceauk/AMR2D_MHD/assets/64028216/d9b8331e-0279-49cc-84d8-463ba2aec858)

3) 2D MHD rotor at t=0.295. 5 levels of grid refinement are used.
   
   a. Density:
   ![rotor045t_0 295062_16x16_lv5_MC_RUSA_density](https://github.com/spaceauk/AMR2D_MHD/assets/64028216/a4022d88-2f0e-473f-8e3f-2b5ce8008e19)
   b. Pressure:
   ![rotor045t_0 295062_16x16_lv5_MC_RUSA_pressure](https://github.com/spaceauk/AMR2D_MHD/assets/64028216/43eedca9-069d-4565-9baf-55835d194585)
   c. Magnetic field along x:
   ![rotor045t_0 295062_16x16_lv5_MC_RUSA_Bx](https://github.com/spaceauk/AMR2D_MHD/assets/64028216/91948d9a-246e-4446-8130-2a8b513a162b)


