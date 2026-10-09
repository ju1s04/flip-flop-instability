PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

COMP       = gnu

USE_MPI    = TRUE

# Native point-mass gravity (castro.use_point_mass) handles gravity; the
# accretion sink is applied through problem_source.H (castro.add_ext_src = 1).
USE_GRAV   = TRUE
USE_REACT  = FALSE

# define the location of the CASTRO top directory
CASTRO_HOME ?= ../../..

# ideal-gas EOS (gamma set at runtime via eos.eos_gamma, matched to problem.gamma_gas)
EOS_DIR     := gamma_law

# no burning -- a single passive species
NETWORK_DIR := general_null
NETWORK_INPUTS = gammalaw.net

PROBLEM_DIR ?= ./

Bpack   := $(PROBLEM_DIR)/Make.package
Blocs   := $(PROBLEM_DIR)

include $(CASTRO_HOME)/Exec/Make.Castro
