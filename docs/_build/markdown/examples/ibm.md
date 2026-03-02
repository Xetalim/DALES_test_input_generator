# Immersed boundary method (IBM) example

The IBM example shows how to represent obstacles (e.g. buildings)
within the domain using the immersed boundary method.

## What this example demonstrates

- Adding an `IBMModule` to a standard atmospheric and LSM setup.
- Using `IBMModification` objects to describe simple geometries such
  as a rectangular building.

## Configuration sketch

The IBM-specific part of the configuration looks like

IBM is typically combined with a grid, atmosphere, LSM and radiation
setup similar to the basic or LSM cases. Once everything is attached,
you run `sim.sim_preprocessing_pipeline()` to generate the IBM input
files.

## When to use this pattern

- Idealised building-resolving experiments.
- Testing how flow responds to simple obstacles before moving to more
  complex geometry import workflows.
  TODO
