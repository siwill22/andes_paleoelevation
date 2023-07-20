41 Reference Elevation Models (REM) are included as sfit objects in BEST_fitting_models.mat. Model elevations can be calculated with each REM in Matlab (Curve Fitting Toolbox is required).

1. Load BEST_fitting_models.mat into Matlab
2. Calculate model elevations h for chemical parameter X at MgO (X and MgO can be scalars or vectors with the same length):
			h = X_model(MgO,X);
	
	e.g., 	h = CeY_CS_model(4.5,5)
			h =
				4.8745
				
			h = CeY_CS_model([2.5 2.8 3.1 4.2],[5 4 3.5 3])
			h = 
				4.5749    4.0321    3.8085    3.7083


	
	