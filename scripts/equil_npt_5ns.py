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
    help="the padding to use in nm",
)
parser.add_argument(
    "-box_size",
    dest="box_size",
    default=None,
    help="the box size to use in nm",
)
parser.add_argument(
    "-n_solvent",
    dest="n_solvent",
    type=int,
    help="the number of solvent particles to use",
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

# Set steps and frequencies
nsteps = 1250000  # 5 ns
report_freq = 100
chk_freq = 500
traj_freq = 2500  # 500 frames

# Setup input
pdb = PDBFile(args.pdb_file)
forcefield = ForceField("amber14/protein.ff14SB.xml", "amber14/tip3p.xml")

# TODO add the dimensions as an option
pdb.topology.setUnitCellDimensions(mm.Vec3(10,10,10) * unit.nanometer)

modeller = Modeller(pdb.topology, pdb.positions)

print("Adding solvent...")

padding = args.padding

# TODO rewrite this part
if not padding:
    print(f" ---> Using custom box vectors and number of solvent atoms")

    if args.box_size is None:
        n_solvent = args.n_solvent
        print(f"number of solvent atoms to use: {n_solvent}")

        modeller.addSolvent(
            forcefield,
            model="tip3p",
            numAdded=n_solvent,
            ionicStrength=0.15 * unit.molar,
        )

    elif args.n_solvent is None:
        box_size = args.box_size

        modeller.addSolvent(
            forcefield,
            model="tip3p",
            boxSize=mm.Vec3(box_size, box_size, box_size),
            ionicStrength=0.15 * unit.molar,
        )

elif padding:
    print(f" ---> Using {geompadding} padding")
    
    geompadding = float(padding) * unit.nanometer

    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=geompadding,
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

# Save the minimised state as a PDB
print(f"Saving minimised state as {state_pdb_filename_min}")
with open(output_prefix + state_pdb_filename_min, "w") as outfile:
    PDBFile.writeFile(
        sim.topology,
        sim.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(),
        file=outfile,
        keepIds=True,
    )

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
