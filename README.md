# CGM_UVB

# UPDATE 26 OCT 2019:

# 1. Grid step size reduced to 0.01.
# 2. Used metallicity = -1 in log scale, i.e, 0.1 in linear scale, instead of 0.1 in log scale.

# ---------------------------------------------------------

# UPDATE 19 OCT 2019:

# 1. File q20_find_nHsol_Zsol.py finds the hydrogen density and metallicity value required by Q20 to match with the observed values, i.e., Q18 data. Changing UVB_Q will lead to values for Q14,Q15, Q16, Q17, Q18, and Q19 too.
# 2. File nHsol_all_elements_plot.py plots out the log(hydrogen density) values by grouping them by the element (C,N,O,S,Si).
# 3. File Zsol_all_elements.py plots out the metallicity values by grouping them by the element (C,N,O,S,Si).

# ---------------------------------------------------------

# UPDATE 17 OCT 2019:

# 2nd file from 13 OCT has been updated to q20_find_nHsol.py which finds the nH at which the column density ratio calculated for Q20 at Z=0.1 matches with the observed column density ratio.
# Results are a bit bizarre and need to be looked into.

# ---------------------------------------------------------

# UPDATE 13 OCT 2019:

# 2 new files added:

# 1. form_obs_col_ratios.py creates a text file storing the column density ratios (14 ratios) of the observations.
# 2. q20_form_colratios.py creates a text file storing column density ratios for Q20, for hydrogen density going from -5 to -3 in 0.05 steps. 

# Next thing to do: Transpose the 41 x 14 array from the second program, and then implement a simple linear interpolation to 
# find nH for each ratio. 

# EXPECTED DATE OF COMPLETION: 19 OCT 2019.

# ---------------------------------------------------------

#
#
#

# q18.py creates a set of test observations using Q18
# and has Z = 0.1 (log scale), hydrogen density = -4 (log scale)
# stopping column density = 10^14 and a redshift of 0.2

# n_H_and_Z_solfinder.py runs an iterative program over a specified range
# of metallicity and hydrogen density for a particular quasar source.
# It lists all possible solution pairs (n_H,Z) for each metal
# In 20 files for 20 metals.
# It also plots these solutions on n_H vs Z space.
