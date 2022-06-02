
TREECOLV2

The TreeCol module is designed to obtain 4pi steradian maps of the sky for EACH gas cell in the simulation. It does this be recording the column density distribution seen by the cell as it walks the gravitational tree. The fact the TreeCol piggy-backs off the tree means that it is locally adaptive in resolution, scales similar to the tree, and comes at only a slight computational cost (around 10% on the treewalk).

The TreeCol results are stored on HealPix maps, which are equal area. The number of pixels in the map is given by NPIX=12*NSIDE*NSIDE, where NSIDE=1, 2, 4, 8, etc. For most purposes, we find that NSIDE=2 is sufficient. 

TreeCol is designed to work with the SGCHEM chemistry module. As such, it can also obtain column densities of H2 and CO, for self-shielding purposes. In principle, it is relatively easy to add more column densities, should they be needed, but see the warning below on memory.

USEAGE

TreeCol is activated using the TREECOLV2 flag in the Config.sh file. The V2 stands for Version 2, since this version of TreeCol uses a lookup table, rather than the method outlined in the original paper (Clark et al. 2012). The lookup table file needs to be present in the run directory of the simulation. The files can be obtained from Paul Clark via email (clarkpc@cardiff.ac.uk).

The main use of TreeCol is to workout the shielding of the interstellar radiation field, and in particular the UV photos responsible for photoelectric emission heating and photo chemistry. Since the photon budget in this energy range is normally dominated by a few nearby stars, we typically want to limit the distance to which we calculate columns densities. We use the parameter TreeColMaxDistance to set the distance out to which TreeCol considers adding particle/node columns to the map. This is code units. 

To obtain H2 and CO maps when using the SGCHEM module, then set flags TREECOLV2_H2 and TREECOLV2_CO respectively. To account for the doppler shifting of the lines in the H2 and CO self-shielding, one can also set the flag TREECOLV2_VEL. This uses the method outlined in Hartwig et al. (2015) to compute the overlap of the lines. TREECOLV2_VEL also requires the parameter FracOverlap to be set in the parameter file, and Hartwig suggests the value 1.694. 

The OUTPUTCOL preprocessor flag will dump the TreeCol maps for EACH CELL in the snapshot (significantly increases the snapshot size!). Note that TreeCol assumes the RING numbering convention for HealPix. This will be important if you ever need to look at the column along a given line of sight. In use with SGCHEM, this is not important, as the ISRF in this module is isotropic.

WARNING: TreeCol significantly increases the memory overhead of the simulation! This is a feature of the Arepo/Gadget family of tree solvers — to find the gravity contribution from particles / cells on other CPUs, the tree sends copies of its local particles to the other MPI tasks, where they then walk the tree belonging to that CPU. This means that all cells/particles (on that timestep) are effectively walking the tree at the same time. We therefore need two copies of the TreeCol maps, one on the home Task, and the other that goes searching for contributions on other tasks. For NSIDE=2, NPIX=48, and so for a basic TreeCol map you need 96 double precision memory addresses PER GAS CELL. For a chemistry run, we need 3 times this, as we store H2 and CO maps as well. So Although TreeCol is relatively CPU cheap, it is memory expensive (it also increases the communication in the tree walk).

PREPROCESSOR FLAGS:
	TREECOLV2
	TREECOLV2_H2    (with SGCHEM)
	TREECOLV2_CO    (with SGCHEM)
	TREECOLV2_VEL   (with SGCHEM and TREECOLV2_H2 and/or TREECOLV2_CO)
	NSIDE=2         (controls the size of the maps and the memory used by TreeCol)
	OUTPUTCOL

Parameters:
	TreeColMaxDistance  (in code units)
	FracOverlap         (with TREECOLV2_VEL)


USEAGE POLICY:

Please contact the authors, Paul Clark (clarkpc@cardiff.ac.uk) and Simon Glover (glover@uni-heidelberg.de) if you wish to use the TreeCol module. In most cases, we will be happy for you to use the module for your publications, without co-authorship, and simply request that you cite the original paper. However, there may be cases where we would request to be involved in the published work.

Original Paper: Clark P.C., Glover S.C.O., Klessen R.S, 2012 MNRAS, 420, 745 (ads link: http://adsabs.harvard.edu/abs/2012MNRAS.420..745C)

