from pathlib import Path

from pymol import cmd

# load structures
pdb_files = sorted(list(Path("../osf/human_abl1").glob("*.pdb")))
for pdb_file in pdb_files:
    cmd.load(pdb_file)

# align structures
for pdb_file in pdb_files[1:]:
    cmd.align(pdb_file.stem, f"{pdb_files[0].stem}")

# protein
cmd.show_as("cartoon", "*")
cmd.color("white", "*")

# highlight kinase motifs
cmd.color("green", "resid 249-254")  # P-loop
cmd.color("red", "resid 282-292")  # aC helix
cmd.color("blue", "resid 381-383")  # DFG
cmd.color("lightblue", "resid 384-402")  # A-loop

# general settings
cmd.bg_colour("white")
cmd.set("ambient", 0.4)
cmd.set("antialias", 1)
cmd.set("ortho", 1)
cmd.set("ray_trace_mode", 1)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_highlight_color", "grey30")

# view
cmd.set_view((
     0.878110528,    0.077458635,   -0.472143412,
     0.275000900,   -0.889225900,    0.365576655,
    -0.391525507,   -0.450857311,   -0.802142560,
    -0.000000075,    0.000114068, -178.545684814,
    27.170295715,    5.755990982,   55.756195068,
   155.570892334,  201.524734497,   20.000000000))

# save pymol session
cmd.save("abl1_structure_superposition.pse")
cmd.ray(5000, 2400)
cmd.png("abl1_structure_superposition.png")
