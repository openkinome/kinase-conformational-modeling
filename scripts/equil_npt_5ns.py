from __future__ import division, print_function

import argparse
import os
import sys
import time
from math import sqrt
from sys import stdout

import mdtraj as md
import openmmtools
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
from mdtraj.reporters import DCDReporter
from simtk.openmm import CustomBondForce, MonteCarloBarostat, XmlSerializer
from simtk.openmm.app import CheckpointReporter, PDBFile
from simtk.openmm.app.forcefield import ForceField
from simtk.openmm.app.modeller import Modeller
from simtk.openmm.app.pdbreporter import PDBReporter
from simtk.openmm.app.statedatareporter import StateDataReporter

parser = argparse.ArgumentParser(description="minimise and equilibrate")
parser.add_argument(
    "-pdb_file",
    dest="pdb_file",
    type=str,
    help="the PDB file to use",
)
parser.add_argument(
    "-padding",
    dest="padding",
    type=float,
    default=0,
    help="the padding to use in nm, default = 0",
)
parser.add_argument(
    "-box_type",
    dest="box_type",
    type=str,
    default="cube",
    help="the box type to use, default = 'cube'",
    choices=["cube", "truncatedOctahedron"],
)
parser.add_argument(
    "-output_dir",
    dest="output_dir",
    type=str,
    default="./",
    help="the output directory to store results",
)

args = parser.parse_args()

# Set parameters
print("Reading parameters...")
pressure = 1.0 * unit.atmospheres
temperature = 310 * unit.kelvin
nonbonded_method = app.PME
constraints = app.HBonds
remove_cm_motion = False
hydrogen_mass = 4.0 * unit.amu  # Using HMR
collision_rate = 1.0 / unit.picoseconds
timestep = 0.004 * unit.picoseconds  # We can use a 4fs timestep with HMR
geompadding = float(args.padding) * unit.nanometer
box_type = args.box_type

# Set steps and frequencies
nsteps = 1250000  # 5 ns
report_freq = 100
chk_freq = 500
traj_freq = 2500  # 500 frames

# Setup input
pdb = PDBFile(args.pdb_file)
forcefield = ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")
modeller = Modeller(pdb.topology, pdb.positions)

print("Adding solvent...")
print(f" ---> Using box type: {box_type}")
print(f" ---> Using {geompadding} padding")

padding, boxSize, boxVectors = None, None, None

if box_type == "cube":

    padding = geompadding

elif box_type == "truncatedOctahedron":

    maxSize = max(
        max((pos[i] for pos in modeller.positions))
        - min((pos[i] for pos in modeller.positions))
        for i in range(3)
    )

    vectors = (
        mm.Vec3(1, 0, 0),
        mm.Vec3(1 / 3, 2 * sqrt(2) / 3, 0),
        mm.Vec3(-1 / 3, sqrt(2) / 3, sqrt(6) / 3),
    )

    boxVectors = [(maxSize + geompadding) * v for v in vectors]


modeller.addSolvent(
    forcefield,
    model="tip3p",
    padding=padding,
    boxVectors=boxVectors,
    boxSize=boxSize,
    ionicStrength=0.15 * unit.molar
)

# Set file names
output_prefix = os.path.join(os.path.normpath(args.output_dir), '') # ensure path string is correct
integrator_xml_filename = "integrator.xml"
state_xml_filename = "state.xml"
state_pdb_filename_min = "min.pdb"
state_pdb_filename = "equilibrated.pdb"
system_xml_filename = "system.xml"
checkpoint_filename = "equilibrated.chk"
traj_output_filename = "equilibrated.dcd"
state_data_filename = "state_data.csv"

# Load the AMBER files
print("Creating OpenMM system input file...")

system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=nonbonded_method,
    constraints=constraints,
    temperature=temperature,
    removeCMMotion=remove_cm_motion,
    hydrogenMass=hydrogen_mass,
)

# Add a barostat to the system
system.addForce(MonteCarloBarostat(pressure, temperature))

# Make and serialize integrator - Langevin dynamics
print(f"Serializing integrator to {integrator_xml_filename}")
integrator = openmmtools.integrators.LangevinIntegrator(
    temperature,
    collision_rate,  # Friction coefficient
    timestep,
    constraint_tolerance=1e-5,
)
with open(output_prefix + integrator_xml_filename, "w") as outfile:
    xml = mm.XmlSerializer.serialize(integrator)
    outfile.write(xml)

# Define the platform to use; CUDA, OpenCL, CPU, or Reference. Or do not specify
# the platform to use the default (fastest) platform
platform = mm.Platform.getPlatformByName("OpenCL")
prop = dict(OpenCLPrecision="mixed")  # Use mixed single/double precision

# Create the Simulation object
sim = app.Simulation(modeller.topology, system, integrator, platform, prop)

# Set the particle positions
sim.context.setPositions(modeller.positions)

print(f"Saving pre-minimised state")
with open(output_prefix + "pre-min.pdb", "w") as outfile:
    PDBFile.writeFile(
        sim.topology,
        sim.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(),
        file=outfile,
        keepIds=True,
    )


# Minimize the energy
print("Minimising energy...")
print(
    "  initial : %8.3f kcal/mol"
    % (
        sim.context.getState(getEnergy=True, enforcePeriodicBox=True).getPotentialEnergy()
        / unit.kilocalories_per_mole
    )
)
sim.minimizeEnergy()
print(
    "  final : %8.3f kcal/mol"
    % (
        sim.context.getState(getEnergy=True, enforcePeriodicBox=True).getPotentialEnergy()
        / unit.kilocalories_per_mole
    )
)

#state = sim.context.getState(getEnergy=True, getForces=True)
#for i, f in enumerate(state.getForces()):
#  if unit.norm(f) > 1e6*unit.kilojoules_per_mole/unit.nanometer:
#    print(i+1, f)

# Save the minimised state as a PDB
print(f"Saving minimised state as {state_pdb_filename_min}")
with open(output_prefix + state_pdb_filename_min, "w") as outfile:
    PDBFile.writeFile(
        sim.topology,
        sim.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(),
        file=outfile,
        keepIds=True,
    )

os.sys.exit()

# set starting velocities:
print("Generating random starting velocities")
sim.context.setVelocitiesToTemperature(temperature)

# write limited state information to standard out:
sim.reporters.append(
    StateDataReporter(
        output_prefix + state_data_filename,
        reportInterval=report_freq,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        temperature=True,
        volume=True,
        speed=True,
        progress=True,
        remainingTime=True,
        totalSteps=nsteps,
        separator="\t",
    )
)

# Write to checkpoint files regularly:
sim.reporters.append(
    CheckpointReporter(
        file=output_prefix + checkpoint_filename, reportInterval=chk_freq
    )
)

# Write out the trajectory
sim.reporters.append(
    md.reporters.DCDReporter(
        file=output_prefix + traj_output_filename, reportInterval=traj_freq
    )
)

# Run NPT dynamics
print("Running dynamics in the NPT ensemble...")
initial_time = time.time()
sim.step(nsteps)
elapsed_time = (time.time() - initial_time) * unit.seconds
simulation_time = nsteps * timestep
print(
    "    Equilibration took %.3f s for %.3f ns (%8.3f ns/day)"
    % (
        elapsed_time / unit.seconds,
        simulation_time / unit.nanoseconds,
        simulation_time / elapsed_time * unit.day / unit.nanoseconds,
    )
)

# Save and serialize the final state
print(f"Serializing state to {state_xml_filename}")
state = sim.context.getState(
    getPositions=True, getVelocities=True, getEnergy=True, getForces=True
)
with open(output_prefix + state_xml_filename, "w") as outfile:
    xml = mm.XmlSerializer.serialize(state)
    outfile.write(xml)

# Save the final state as a PDB
print(f"Saving final state as {state_pdb_filename}")
with open(output_prefix + state_pdb_filename, "w") as outfile:
    PDBFile.writeFile(
        sim.topology,
        sim.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(),
        file=outfile,
        keepIds=True,
    )

# Save and serialize system
print(f"Serializing system to {system_xml_filename}")
system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
with open(output_prefix + system_xml_filename, "w") as outfile:
    xml = mm.XmlSerializer.serialize(system)
    outfile.write(xml)
